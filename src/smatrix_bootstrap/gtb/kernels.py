"""Kernel export format (shared by the Mathematica and Arb generators): writer and loader.

An export is a directory with ``manifest.json`` and one ``<name>.txt.gz`` per matrix:
row-major, one decimal per line written as ``<int>e<exp>`` (``digits`` significant
digits) or ``0`` for an exact zero.  The manifest records
``{"M", "L", "nu0", "digits", "precision", "generator", "matrices": {name: {"shape",
"file", "sha256"}}}``; ``sha256`` is the digest of the stored ``.gz`` bytes (the digest
of the decompressed text is also accepted, so either generator convention loads).
"""
from __future__ import annotations

import gzip
import hashlib
import json
from math import prod
from pathlib import Path

import numpy as np


class KernelSet:
    """A verified export: ``ks.<name>`` is a float64 array, ``ks.exact(name)`` the decimals."""

    def __init__(self, path: Path, manifest: dict):
        self.path = path
        self.manifest = manifest
        self._cache: dict[str, np.ndarray] = {}

    @property
    def names(self) -> list[str]:
        return list(self.manifest["matrices"])

    def _text(self, name: str) -> str:
        meta = self.manifest["matrices"][name]
        return gzip.decompress((self.path / meta["file"]).read_bytes()).decode()

    def exact(self, name: str) -> list[str]:
        lines = self._text(name).split()
        if len(lines) != prod(self.manifest["matrices"][name]["shape"]):
            raise ValueError(f"{name}: entry count does not match the manifest shape")
        return lines

    def array(self, name: str) -> np.ndarray:
        if name not in self._cache:
            shape = tuple(self.manifest["matrices"][name]["shape"])
            text = self._text(name)
            arr = np.array(text.split(), dtype=np.float64)
            if arr.size != prod(shape):
                raise ValueError(f"{name}: {arr.size} entries, manifest shape {shape}")
            self._cache[name] = arr.reshape(shape)
        return self._cache[name]

    def __getattr__(self, name: str):
        if name.startswith("_") or name not in self.manifest.get("matrices", {}):
            raise AttributeError(name)
        return self.array(name)


def load(path) -> KernelSet:
    """Read ``manifest.json`` and check every file's sha256 (arrays are parsed lazily)."""
    path = Path(path)
    manifest = json.loads((path / "manifest.json").read_text())
    for name, meta in manifest["matrices"].items():
        raw = (path / meta["file"]).read_bytes()
        digests = {hashlib.sha256(raw).hexdigest()}
        if meta["sha256"] not in digests:
            digests.add(hashlib.sha256(gzip.decompress(raw)).hexdigest())
        if meta["sha256"] not in digests:
            raise ValueError(f"{name}: sha256 mismatch for {meta['file']}")
    return KernelSet(path, manifest)


def fmt(v, digits: int = 40, stats: dict | None = None) -> str:
    """An Arb ball as ``<int>e<exp>`` (midpoint, ``digits`` significant), ``0`` if it holds 0."""
    if 0 in v:
        if stats is not None and not v.is_zero():
            stats["zero_balls"] = stats.get("zero_balls", 0) + 1
        return "0"
    mant, _, ex = v.mid().str(digits, radius=False).partition("e")
    ip, _, fp = mant.partition(".")
    sign, ip = ("-", ip[1:]) if ip.startswith("-") else ("", ip)
    body = (ip + fp).lstrip("0")
    digits_ = body.rstrip("0")
    return f"{sign}{digits_}e{int(ex or 0) - len(fp) + len(body) - len(digits_)}"


class _Handle:
    def __init__(self, writer: "ExportWriter", name: str, shape):
        self.writer, self.name, self.shape = writer, name, [int(n) for n in shape]
        self.file = writer.outdir / f"{name}.txt.gz"
        self._raw = open(self.file, "wb")
        self._gz = gzip.GzipFile(filename="", mode="wb", fileobj=self._raw, mtime=0)
        self.count = 0

    def add(self, lines) -> None:
        if lines:
            self._gz.write(("\n".join(lines) + "\n").encode())
            self.count += len(lines)

    def close(self) -> None:
        self._gz.close()
        self._raw.close()
        if self.count != prod(self.shape):
            self.file.unlink()
            raise ValueError(f"{self.name}: wrote {self.count} entries, shape {self.shape}")
        digest = hashlib.sha256(self.file.read_bytes()).hexdigest()
        self.writer.entries[self.name] = {"shape": self.shape, "file": self.file.name,
                                          "sha256": digest}


class ExportWriter:
    """Streamed writer: ``open(name, shape).add(row_strings)...close()``, then ``finish()``."""

    def __init__(self, outdir, meta: dict):
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.meta = dict(meta)
        self.entries: dict[str, dict] = {}

    def open(self, name: str, shape) -> _Handle:
        return _Handle(self, name, shape)

    def write(self, name: str, shape, rows) -> None:
        h = self.open(name, shape)
        for row in rows:
            h.add(list(row))
        h.close()

    def finish(self, **extra) -> dict:
        manifest = {**self.meta, **extra, "matrices": dict(sorted(self.entries.items()))}
        (self.outdir / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
        return manifest
