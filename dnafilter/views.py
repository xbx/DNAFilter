from django.shortcuts import render_to_response, get_object_or_404
from django.template import RequestContext
from dnafilter.models import SequenceBase

from Bio import SeqIO
from filtro3g import FiltroSec, PaintSeq

import sys
sys.stdout = sys.stderr

def is_fasta(text):
    for count, line in enumerate(text.split('\n')):
        if not line:
            return False
        if count == 0 and line[0] != '>':
            return False
        if line[0] != '>' and ' ' in line.strip():
            return False
    return True


def index(request):
    db_list = SequenceBase.objects.all()
    # TO DO:
    # Process db_list to remove path (path disclosure).
    return render_to_response('index.html', {'db_list': db_list}, 
                              context_instance=RequestContext(request))

def filter(request):
    
    # TO DO: READ FROM FILE
    cfg = {'paths': 
                     {
                      'makeblastdb_exe': 'makeblastdb', 
                      'blast_exe': 'blastn', 
                      'db_path': '/var/www/dnafilter/db'
                     }
           }
    
    
    form = {
        'blastdb':      request.POST.get('filter-by'),
        'seqs':         request.POST.get('seqs'),
        'seqdatafile':  request.POST.get("seqdatafile"),
        'input_method': request.POST.get('input_method'),
        'order':        request.POST.get("order"),
        'connector':    request.POST.get("connector"),
        'linker1':      request.POST.get("linker1"),
        'linker2':      request.POST.get("linker2"),
        'isize':        request.POST.get("isize", '5'),
        'colors': {
            'vector':       request.POST.get("vector-color"),
            'connector':    request.POST.get("connector-color"),
            'linker':       request.POST.get("linker-color"),
            'mirna':          request.POST.get("mirna-color"),
        }
    }
    
    # determine if we're using direct input or uploaded file
    if form['input_method'] == 'seqdatafile':
        seqs = ''
        if form['seqdatafile'] is not None:
            seqs_f = form['seqdatafile'].file
            try:
                #TODO handle errors (LargeZipFile e, z.testzip)
                z = ZipFile(seqs_f)
                seqs = '\n'.join(z.open(z_fname).read() for z_fname in z.namelist() if not z_fname.endswith('/'))
            except BadZipfile:
                seqs_f.seek(0)
                seqs = seqs_f.read()
    else:
        seqs = form['seqs']
        if not is_fasta(seqs):
            return render_to_response('index.html', {'error': True})
            
    
    # create instance of filter-tool and put uploaded seqs in a temporary file
    filtro = FiltroSec(cfg)
    fasta_in_fn = filtro.filesfromseqs(seqs)

    # apply filters with the specified database
    filtered_fn, bat1_col = filtro.apply_filter(fasta_in_fn, form['blastdb'], 
                                                color='vector')
                                                
    # create instance of filter-tool and put uploaded seqs in a temporary file
    filtro = FiltroSec(cfg)
    fasta_in_fn = filtro.filesfromseqs(seqs)

    # apply filters with the specified database
    filtered_fn, bat1_col = filtro.apply_filter(fasta_in_fn, form['blastdb'], 
                                                color='vector')

    # Do not color if bat1_col is empty
    col_seqs = []
    if bat1_col:
        for seq in SeqIO.parse(fasta_in_fn, 'fasta'):
            col_seq = PaintSeq(seq.seq, bat1_col[seq.id], 50).new_seq_html
            col_seqs.append('>%s\n%s' % (seq.id, col_seq))
    else:
        for seq in SeqIO.parse(fasta_in_fn, 'fasta'):
            col_seqs.append('>%s\n%s' % (seq.id, seq.seq))
    form['finalout_color'] = '<br>'.join(col_seqs)
    form['finalout'] = open(filtered_fn).read()

    fasta_cut_fn = filtro.cut(filtered_fn)
    with open(fasta_cut_fn) as finalcut:
        form['finalcut'] = finalcut.read()


    print 'fasta_cut_fn: %s' % form['finalcut']

    if form['connector']:
        # Filtered sequence: fasta_cut_fn
        connector_seq = '>Conector\n%s' % form['connector']
        connector_seq_fn = filtro.filesfromseqs(connector_seq)
        connector_fil_fn, bat1_col = filtro.apply_filter(fasta_in_fn,
                                               connector_seq_fn,
                                               'connector', 'O', 
                                               'connector', bat1_col)
        
    if form['linker1'] and form['linker2']:
        lnk1_seq = ">Linker 3'\n%s" % form['linker1']
        lnk1_seq_fn = filtro.filesfromseqs(lnk1_seq)
        lnk2_seq = ">Linker 5'\n%s" % form['linker2']
        lnk2_seq_fn = filtro.filesfromseqs(lnk2_seq)

        lnk1_fil_fn, bat1_col = filtro.apply_filter(connector_fil_fn,
                                          lnk1_seq_fn, 'lnk', 'L',
                                          'linker', bat1_col)

        lnk2_fil_fn, bat1_col = filtro.apply_filter(lnk1_fil_fn,
                                          lnk2_seq_fn, 'lnk','L',
                                          'linker', bat1_col)

    col_seqs = []
    targets_seqs = []
    i = 0
    for seq in SeqIO.parse(fasta_in_fn, 'fasta'):
        i += 1
        col_seq = PaintSeq(seq.seq, bat1_col[seq.description], 50, form['isize'])
        col_seqs.append('>%s\n%s' % (seq.id, col_seq.new_seq_html))
        targets = []
        j = 0
        for target_seqs_idx in col_seq.target_seqs_idxs:
            j += 1
            target_seq = str(seq.seq[target_seqs_idx[0]:
                                     target_seqs_idx[1]+1])
            targets.append('>Insert#%s from %s\n%s' % (j, seq.id,
                                                      target_seq))
        targets_html = '<br>'.join(targets)
        targets_seqs.append(targets_html)
        #col_seq.target_seqs_idxs
        #targets_seqs.append(col_seq.target_seqs_idxs)
        #print dir(seq), '*************************'
    form['finalout_color2'] = '<br>'.join(col_seqs)
    form['target_seqs'] = '<br>'.join(targets_seqs)
    
    
    return render_to_response('filter.html', {'form': form})
#                              context_instance=RequestContext(request))
