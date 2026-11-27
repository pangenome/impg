#!/usr/bin/env python3
"""Public five-assembly AGC + real SYNG + native junction-query integration gate. No specimen data."""
import argparse
import collections
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess as sp
import tempfile

RC=str.maketrans('ACGT','TGCA')
def rc(s): return s.translate(RC)[::-1]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x): p.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')

def fixture(impg, agc, root):
    root.mkdir(exist_ok=False); rng=random.Random(59631)
    seq=lambda n: ''.join(rng.choices('ACGT',k=n))
    a,b=seq(36000),seq(34100); p,q=5400,9600; x,y=15400,17600
    # Two transferable bases in a second (zero-insert) model.
    b=b[:y]+a[x:x+2]+b[y+2:]; rep,rep2=21000,22000
    a=a[:30000]+'AT'*110+a[30220:];a=a[:22700]+'N'*800+a[23500:]
    c1='REF#0#one';c2='REF#0#two'; lib={}
    for g in ('REF','B','C','D'):
        # Both arms exist in this assembly, but on two contigs, not joined.
        lib[g+'#0']={g+'#0#one':a,g+'#0#two':b}
    e=(seq(13000)+a[p-200:p]+'AGC'+b[q:q+200]+seq(450)+rc(a[x-200:x]+b[y:y+200])+seq(450)
       +a[rep-200:rep]+b[rep2:rep2+200]+seq(450)+a[rep-200:rep]+seq(15000))
    lib['E#1']={'E#1#witness':e}
    def m(id,c,k,s,d,h,t,ins,eq=None):return dict(id=id,end1=[c,k,s],end2=[d,h,t],insert=ins,coupled=eq or [[k,h]])
    ms=[m('fifth',c1,p,'L',c2,q,'R','AGC'),m('reverse',c2,q,'R',c1,p,'L','GCT'),
        m('micro',c1,x,'L',c2,y,'R','',[[x,y],[x+1,y+1],[x+2,y+2]]),
        m('repeat',c1,rep,'L',c2,rep2,'R',''),m('wrong_orientation',c1,p,'L',c2,q,'L','AGC'),
        m('wrong_insert',c1,p,'L',c2,q,'R','AAC'),m('unknown_insert',c1,p,'L',c2,q,'R',None),
        m('edge',c1,3,'L',c2,q,'R',''),m('masked',c1,22900,'L',c2,q,'R',''),
        m('overlap_palindrome',c1,30120,'L',c1,30120,'R','')]
    combined=root/'panel.fa'; fas=[]
    with combined.open('w') as out:
        for group,paths in lib.items():
            fp=root/(group.split('#')[0]+'.fa');fas.append(fp)
            with fp.open('w') as f:
                for name,s in paths.items():
                    record='>'+name+'\n'+s+'\n';f.write(record);out.write(record)
    archive=root/'panel.agc';prefix=root/'graph';ledger=root/'models.jsonl';ledger.write_text(''.join(json.dumps(z)+'\n' for z in ms))
    sp.run([agc,'create','-t','2','-o',archive,*fas],check=True)
    sp.run([impg,'syng','-t','2','-v','0','--syncmer-length','63','--smer-length','8','--syncmer-seed','7','--position-sample-rate','64','-f',combined,'-o',prefix],check=True)
    def run(out,extra=(),check=True,models=None):
        return sp.run([impg,'junction-query','--index',prefix,'--assemblies',archive,'--models',models or ledger,'--out',out,'-t','2','--examples','4',*extra],check=check)
    roles={'flank1_40':(40,1),'flank2_40':(40,2),'join22':(22,0),'join33':(33,0),'joined40':(40,0),'flank1_96':(96,1),'flank2_96':(96,2),'joined96':(96,0)}
    def literal(mod,w,side):
        def arm(end,left):
            c,k,t=end;ref=lib['REF#0'][c];z=ref[max(0,k-w):k] if t=='L' else ref[k:k+w];return z if (t=='L')==left else rc(z)
        left=arm(mod['end1'],True);right=arm(mod['end2'],False)
        if side==1:v=left if len(left)==w else None
        elif side==2:v=right if len(right)==w else None
        else:v=left+mod['insert']+right if mod['insert'] is not None and len(left)==len(right)==w else None
        return v if v is not None and set(v)<=set('ACGT') else None
    def count(sequence,pattern):
        n=0;i=sequence.find(pattern)
        while i>=0:n+=1;i=sequence.find(pattern,i+1)
        return n
    for name,sets in [('four',list(lib)[:-1]),('five',list(lib))]:
        out=root/name;run(out,['--sets',','.join(sets)]);summ=[json.loads(l) for l in (out/'summary.jsonl').read_text().splitlines()];assert [z['id'] for z in summ]==[z['id'] for z in ms]
        for z,mod in zip(summ,ms):
            for k,(w,side) in roles.items():
                v=z['literal'][k];word=literal(mod,w,side)
                if word is None:assert v['status']=='not_assessed' and v['occurrences'] is None;continue
                gold={g:sum(count(s,t) for s in lib[g].values() for t in {word,rc(word)}) for g in sets};gold={g:n for g,n in gold.items() if n}
                assert v['by_set']==gold and v['occurrences']==sum(gold.values()) and v['status']==('present' if gold else 'absent_after_exhaustive_search'),(k,mod['id'])
        by={z['id']:z for z in summ};assert by['fifth']['placement']==('verified_same_path_join' if name=='five' else 'no_verified_join');assert by['wrong_insert']['placement']==by['wrong_orientation']['placement']=='no_verified_join';assert by['unknown_insert']['placement']=='unassessable'
        if name=='five':
            assert by['repeat']['placement']=='ambiguous_locus' and by['micro']['unique_anchor_sets']==['E#1'];hits=[json.loads(l) for l in (out/'hits.jsonl').read_text().splitlines()];f,=[r for r in hits if r['id']=='fifth'];bck,=[r for r in hits if r['id']=='reverse'];mh,=[r for r in hits if r['id']=='micro'];assert f['strand']=='+' and bck['strand']=='-' and f['arm1']==bck['arm2'] and f['arm2']==bck['arm1'];assert mh['strand']=='-' and mh['coupled_target_cuts0'][2]==[x-2 for x in mh['coupled_target_cuts0'][0]]
        c=json.loads((out/'complete.json').read_bytes());assert c['paths_assessed']==sum(len(lib[g]) for g in sets) and c['bases_assessed']==sum(len(s) for g in sets for s in lib[g].values()) and c['signed_midpath_steps_compared_to_source']>0 and c['summary_sha256']==sha(out/'summary.jsonl')
    out=root/'five';before=sha(out/'summary.jsonl');run(out,['--sets',','.join(lib)]);assert sha(out/'summary.jsonl')==before
    bad=root/'bad.jsonl';bad.write_text(json.dumps(dict(ms[2],coupled=[[x,y],[x+1,y]]))+'\n');assert run(root/'invalid',models=bad,check=False).returncode!=0 and not (root/'invalid/complete.json').exists()
    unit=next((out/'units').glob('*.json'));unit.write_bytes(unit.read_bytes()+b' ');assert run(out,['--sets',','.join(lib)],check=False).returncode!=0;unit.write_bytes(unit.read_bytes()[:-1]);run(out,['--sets',','.join(lib)])
    result={'case':'public synthetic; exact per-role/set naive overlap oracle; held-out fifth path, both orientations, wrong insert/strand, same hap different contigs, repeat, N, edge, unknown insert, palindrome, coupled shift and corrupt seal','passed':True,'model_count':len(ms),'input_sha256':sha(ledger),'AGC_sha256':sha(archive),'tool_sha256':sha(Path(impg)),'outputs':{k:sha(root/k/'complete.json') for k in ('four','five')}};dump(root/'fixture.json',result);return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('impg');p.add_argument('agc');p.add_argument('output',type=Path);args=p.parse_args();fixture(args.impg,args.agc,args.output);print('PUBLIC_NATIVE_JUNCTION_FIXTURE_PASS',flush=True)
