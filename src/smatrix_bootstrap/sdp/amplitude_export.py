"""Exact data conversion of completed SDPB coordinates to density amplitudes."""
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import math
import re
import time

import numpy as np


def _digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def _decimal_fraction(text):
    if len(text) > 4096 or not re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?',text):
        raise ValueError('Expected a bounded finite decimal SDPB token')
    try:
        value = Decimal(text)
    except InvalidOperation as exc:
        raise ValueError('Invalid decimal SDPB token') from exc
    if not value.is_finite() or abs(value.as_tuple().exponent) > 10000:
        raise ValueError('Decimal exponent outside the exact-export limit')
    return Fraction(value)


def exact_coefficients(basis, a_text, *, deadline=None, max_entries=50_000_000):
    """Exact integer dot products, using one dyadic denominator per basis row."""
    a = [_decimal_fraction(value) for value in a_text]
    if basis is None:
        return a
    if (basis.ndim != 2 or basis.shape[1] != len(a) or basis.size > max_entries or
            basis.dtype != np.dtype('float64') or not np.all(np.isfinite(basis))):
        raise ValueError('Expected a finite binary64 basis within the entry limit')
    denominator = math.lcm(*(value.denominator for value in a))
    integers = [value.numerator*(denominator//value.denominator) for value in a]
    out = []
    for row in basis:
        if deadline is not None and time.monotonic() > deadline:
            raise TimeoutError('Exact export deadline reached between basis rows')
        entries = [float(value).as_integer_ratio() for value in row]
        binary_denominator = max(q for _,q in entries)
        numerator = sum(p*(binary_denominator//q)*v for (p,q),v in zip(entries,integers))
        out.append(Fraction(numerator,binary_denominator*denominator))
    return out


def _source(source):
    """Authenticate finite output data without evaluating or assembling an SDP."""
    from .projector import Layout
    rec = json.loads(source.read_text());root = source.parent
    if (rec.get('solver') != 'SDPB' or rec.get('status') not in ('numerically_accepted','not_accepted')
            or rec.get('sdpb',{}).get('returncode') != 0 or 'verification' not in rec
            or any(key in rec for key in ('solver_failure','readback_error','untrusted_partial_output'))):
        raise ValueError('Export requires a completed SDPB leaf with zero solver exit status')
    process_path = root/'sdpb_process.json'
    process = json.loads(process_path.read_text())
    if process.get('terminal') is not True or process.get('returncode') != 0:
        raise ValueError('Solver process is unfinished or exited unsuccessfully')
    reason = rec.get('sdpb_result',{}).get('terminateReason','')
    if (not reason or re.search(r'time|runtime|max|cancel|interrupt|signal|terminated',reason,re.I) or
            rec.get('convergence',{}).get('terminate_reason') != reason):
        raise ValueError('Missing/conflicting terminal reason or budget/termination residual output')
    spec,pmp = rec['spec'],rec['pmp'];M,L,na = spec['M'],spec['L'],pmp['n_a']
    if any(type(value) is not int or value <= 0 for value in (M,L,na)):
        raise ValueError('Invalid amplitude dimensions')
    layout = Layout(M);ny = na+(4*M if spec['uv'] else 0)
    if (pmp['n_vars'] != ny+1 or na > layout.n or pmp.get('reduce_basis') != spec['reduce_basis']
            or (not spec['reduce_basis'] and na != layout.n)):
        raise ValueError('Recorded coordinate layout mismatch')
    scheme = spec.get('scattering_prescription','mixed-pv')
    if pmp.get('operator_precision',{}).get('scattering_prescription',scheme) != scheme:
        raise ValueError('Recorded source interpolation identity mismatch')
    for name in ('projector.py','grid.py'):
        if rec.get('source_sha256',{}).get('src/smatrix_bootstrap/sdp/'+name) != _digest(Path(__file__).with_name(name)):
            raise ValueError('Source density layout/node producer changed: '+name)
    files = [source,root/'pmp.json',root/'out/y.txt',root/'solution.npz',process_path]
    basis = None
    if spec['reduce_basis']:
        files.append(root/'basis.npy')
        basis = np.load(files[-1],mmap_mode='r',allow_pickle=False)
        if basis.shape != (layout.n,na) or basis.dtype != np.dtype('float64'):
            raise ValueError('Saved basis dimensions/type differ from the source layout')
    hashes = {str(path):_digest(path) for path in files}
    expected = {'pmp.json':pmp['sha256'],'solution.npz':rec.get('solution_sha256'),
                'out/y.txt':rec.get('y_sha256')}
    if basis is not None:
        expected['basis.npy'] = rec.get('basis',{}).get('sha256')
    for name,digest in expected.items():
        if digest is None or digest != hashes[str(root/name)]:
            raise ValueError('Missing solve-time hash or source input hash mismatch: '+name)
    lines = (root/'out/y.txt').read_text().splitlines()
    if not lines or lines[0].split() != [str(ny),'1']:
        raise ValueError('Full y text header differs from the coordinate layout')
    tokens = [line.strip() for line in lines[1:] if line.strip()]
    if len(tokens) != ny:
        raise ValueError('Incomplete full y text')
    for token in tokens:
        _decimal_fraction(token)
    rounded = np.array([float(token) for token in tokens])
    if not np.all(np.isfinite(rounded)):
        raise ValueError('Saved y cannot match a finite binary64 solution copy')
    with np.load(root/'solution.npz',allow_pickle=False) as saved:
        arrays = {'y':rounded,'a':rounded[:na]}
        if spec['uv']:
            arrays.update(ImF=rounded[na:na+2*M].reshape(2,M),rho_hat=rounded[na+2*M:].reshape(2,M))
        if saved['c'].shape != (layout.n,) or any(key not in saved or not np.array_equal(saved[key],value)
                                                for key,value in arrays.items()):
            raise ValueError('Authenticated NPZ shape/y/current copy mismatch')
    return rec,layout,basis,tokens,hashes


def _blocks(c, layout):
    M = layout.M
    return {'T0':c[0],'sigma1':c[layout.s1],'sigma2':c[layout.s2],
        'rho1_row_major':[c[layout.r1][i*M:(i+1)*M] for i in range(M)],
        'rho2_upper_unscaled':c[layout.r2]}


def _strings(value):
    if isinstance(value,Fraction):return str(value)
    if isinstance(value,dict):return {key:_strings(item) for key,item in value.items()}
    if isinstance(value,list):return [_strings(item) for item in value]
    return value


def _wolfram(value):
    """Only literals, Lists and Associations: no calls, assignments or side effects."""
    if isinstance(value,Fraction):return str(value)
    if isinstance(value,int):return str(value)
    if isinstance(value,str):return json.dumps(value,ensure_ascii=True)
    if isinstance(value,list):return '{'+','.join(_wolfram(item) for item in value)+'}'
    if isinstance(value,dict):return '<|'+','.join(_wolfram(key)+'->'+_wolfram(item) for key,item in value.items())+'|>'
    raise TypeError('Nonliteral Wolfram export value')


def export_amplitude(source_report,outdir,timeout=600,max_entries=50_000_000):
    """Export the exact amplitude represented by authenticated saved finite data."""
    started = time.monotonic();source,dest = Path(source_report).resolve(),Path(outdir).resolve()
    if dest.exists() and (not dest.is_dir() or any(dest.iterdir())):
        raise FileExistsError('Choose a fresh exact-amplitude export directory')
    if not math.isfinite(timeout) or timeout <= 0 or max_entries <= 0:
        raise ValueError('Positive finite timeout and basis entry limit required')
    rec,layout,basis,tokens,hashes = _source(source);na = rec['pmp']['n_a'];M = layout.M
    c = exact_coefficients(basis,tokens[:na],deadline=started+timeout,max_entries=max_entries)
    blocks = _blocks(c,layout)
    pairs = [[int(i),int(j)] for i,j in zip(*layout.triu)]
    currents = None
    if rec['spec']['uv']:
        currents = {'ImF':{wave:tokens[na+ell*M:na+(ell+1)*M] for ell,wave in enumerate(('S0','P1'))},
            'rho_hat':{wave:tokens[na+2*M+ell*M:na+2*M+(ell+1)*M] for ell,wave in enumerate(('S0','P1'))},
            'representation':'original full y decimal token strings, unchanged',
            'rho_definition':'rho_ell(s)=k_ell(s)^2 rho_hat_ell(s)',
            'k_squared':{'S0':'3 sqrt(1-4/s)/(256 pi^5)','P1':'(s-4) sqrt(1-4/s)/(384 pi^5)'},
            'form_factor':'ImF is stored; ReF_i=1+sum_j K_ij ImF_j; F(0)=1',
            'K_definition':'K_ij=-v[(i-j) mod 2M]+v[(i+j+1) mod 2M], zero-based; v[n]=cot(pi n/(2M))/M for odd n, otherwise 0'}
    payload = {'schema':'canonical_density_amplitude_v1','source_report':str(source),
        'source_sha256':hashes[str(source)],'input_sha256':hashes,'spec':rec['spec'],
        'recorded_source_sha256':rec.get('source_sha256',{}),'source_numerically_accepted':rec.get('accepted',False),
        'solver_terminal_reason':rec['sdpb_result']['terminateReason'],
        'representation':'reduced exact rational strings p/q or integer; every value is exact for the saved finite inputs',
        'conversion':{'method':'full decimal y -> Fraction; float.as_integer_ratio basis -> dyadics; common-denominator integer dot product per row; final exact reduction',
            'representation_error':'0 exactly','basis_memory':'binary64 memory map; one row of integer ratios at a time',
            'basis_entries':int(basis.size) if basis is not None else 0,
            'max_basis_entries':max_entries,'max_decimal_token_characters':4096,'max_absolute_decimal_exponent':10000,
            'timeout_seconds':timeout,'deadline_scope':'checked between basis rows'},
        'layout':{'M':M,'L':rec['spec']['L'],'n_c':layout.n,'n_a':na,'n_y':len(tokens),
            'source_y_order':'a, ImF_S0, ImF_P1, rho_hat_S0, rho_hat_P1; normalized y0 omitted; current blocks absent for IR',
            'c_order':'T0, sigma1, sigma2, rho1 row-major, rho2 upper triangle row-major including diagonal',
            'rho2_upper_pairs_zero_based':pairs,'rho2_unpack':'rho2[i,j]=rho2[j,i]=packed value; no factor of 2 or sqrt(2)'},
        'node_convention':{'index':'j=0,...,M-1','phi_j':'pi (2j+1)/(2M)',
            's_j':'8/(1+cos(phi_j)); dimensionless m_pi=1','energy_GeV':'(140/1000) sqrt(s_j)',
            'scattering_prescription':rec['spec'].get('scattering_prescription','mixed-pv')},
        'amplitude':_strings(blocks),'currents':currents,
        'certified':False,'optimization_performed':False,
        'scope':'exact representation conversion of one saved finite amplitude; SVD/model truncation, solver accuracy, global unitarity and physics reproduction are not certified or changed',
        'producer_sha256':{name:_digest(Path(__file__).with_name(name)) for name in ('amplitude_export.py','projector.py','grid.py')}}
    rho2 = [[Fraction(0) for _ in range(M)] for _ in range(M)]
    for (i,j),value in zip(pairs,blocks['rho2_upper_unscaled']):rho2[i][j]=rho2[j][i]=value
    wl = {'schema':'canonical_density_amplitude_v1','M':M,'L':rec['spec']['L'],
          'scattering_prescription':payload['node_convention']['scattering_prescription'],
          **blocks,'rho2_full_symmetric':rho2}
    if currents is not None:
        for key in ('ImF','rho_hat'):
            wl[key] = {wave:[_decimal_fraction(token) for token in values] for wave,values in currents[key].items()}
    if any(_digest(path) != digest for path,digest in hashes.items()):
        raise ValueError('An authenticated source file changed during export')
    payload['conversion']['seconds'] = time.monotonic()-started
    dest.mkdir(parents=True,exist_ok=True)
    if any(dest.iterdir()):raise FileExistsError('Export directory became nonempty')
    (dest/'amplitude.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    (dest/'amplitude.wl').write_text(_wolfram(wl)+'\n')
    (dest/'canonical_y.txt').write_bytes((source.parent/'out/y.txt').read_bytes())
    report = {'kind':'canonical_amplitude_export','status':'complete','source_report':str(source),
        'source_sha256':hashes[str(source)],'n_c':layout.n,'n_y':len(tokens),
        'seconds':time.monotonic()-started,'representation_error':'0 exactly',
        'files_sha256':{name:_digest(dest/name) for name in ('amplitude.json','amplitude.wl','canonical_y.txt')},
        'producer_sha256':payload['producer_sha256'],'scope':payload['scope'],'certified':False}
    (dest/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
