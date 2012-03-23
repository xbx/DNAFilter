#!/usr/bin/env python

# Revision: $Rev:$

from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline as blastCL

from tempfile import NamedTemporaryFile
import copy
import sys
import os
import json
import itertools
import time

def isoverlapped(s1, s2):

    def range_order(x):
        if x[0]>x[1]:
            return (x[1], x[0])
        else:
            return x

    s1 = range_order(s1)
    s2 = range_order(s2)

    if s1[0]<=s2[0]<=s2[1]<=s1[1]:
        return True
    elif s1[0]<=s2[0]<=s1[1]<=s2[1]:
        return True
    elif s2[0]<=s1[0]<=s2[1]<=s1[1]:
        return True
    elif s2[0]<=s1[0]<=s1[1]<=s2[1]:
        return True
    return False

def extended_range(s1, s2):
    if s1[0]<=s2[0]<=s2[1]<=s1[1]:
        return s1
    elif s2[0]<=s1[0]<=s1[1]<=s2[1]:
        return s2
    elif s1[0]<=s2[0]<=s1[1]<=s2[1]:
        return (s1[0], s2[1])
    elif s2[0]<=s1[0]<=s2[1]<=s1[1]:
        return (s2[0], s1[1])


class FiltroSec(object):
    """

    Test:

    >>> 2 + 2
    4

    >>> import json
    >>> f = FiltroSec(json.load(open("config.json")))
    >>> fhf = f.apply_filter()
    >>> x = open(fhf).read()
    >>> x == f._testseqOK
    False
    """

    _testseq = '''>090819-17_E04_129-M13F.ab1   950
NNNNNNNNNNTAGGGCGATTGATTTAGCGGCCGCGAATTCGCCCTTGATT
GATGGTGCCTACAGGTGTGTACAAAGGGCAGGGACCTTGGTGCCTACAGA
GACATTCGAAACACTACCTCAACCTTGGTGCCTACAGGGATACCTTGGTG
CCCGAGAATTCCAAAGGGCGAATTCGTTTAAACCTGCAGGACTAGTCCCT
TTAGTGAGGGTTAATTCTGAGCTTGGCGTAATCATGGTCATAGCTGTTTC
CTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGA
AGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATT
AATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCC
AGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATT
GGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCG
GCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCA
CAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAA
AAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCT
CCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGC
GAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCC
CTCGTGCGCTCTGCTGTTCCGACCCTGCCGCTTTACCGGGATACCTGTCC
GCCTTTCTCCCTTCGGGAAGCGTGCGCTTTCTCATATCTCACGCTGTAGG
TATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACAA
CCCCCCGGTCAGCCCGCCGCTGCGCCTTATCCGGTATGTATCGACTTTGG
'''

    _testseqOK = '''>090819-17_E04_129-M13F.ab1 950
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
GATGGTGCCTACAGGTGTGTACAAAGGGCAGGGACCTTGGTGCCTACAGA
GACATTCGAAACACTACCTCAACCTTGGTGCCTACAGGGATACCTTGGTG
CCCGAGAATTCCAAAGGGCGAATTCGTTTAAACCTGCAGGACTAGTCCCT
TTAGTGAGGGTTAATTCTGAGCTTGGCGTAATCATGGTCATAGCTGTTTC
CTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGA
AGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATT
AATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCC
AGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATT
GGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCG
GCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCA
CAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAA
AAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCT
CCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGC
GAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCC
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
'''



    def __init__(self, cfg):
        """
        Permite indicar las rutas de acceso a la base de vectores
        utilizada para la busqueda de coincidencias y
        la herramienta BLAST para realizar el alineamiento
        TODO: Param. mask character.
        """
        self.blastn_path = cfg['paths']["blast_exe"]
        self.__formato_input = "fasta"
        self.__mask = "N"

    def filesfromseqs(self, seq):
        """ Save a seq object into a file with a tmp name.
        """
        # TODO: Check if it is in Seq format.
        with NamedTemporaryFile(delete=False) as fasta_in_fh:
            fasta_in_fh.write(seq)
            return fasta_in_fh.name

    def sacNfinales(self, nroN, seqI):
        # desde la derecha.
        if 'N' in seqI[-nroN:]:
            # busco el de mas a la izq:
            cola = seqI[-nroN:]
            npos = cola.index('N')
            seqI = seqI[:-nroN+npos]
            seqI = self.sacNfinales(nroN, seqI)
        return seqI

    def sacNiniciales(self, nroN, seqI):
        # desde la izq.
        if 'N' in seqI[:nroN]:
            # busco el de mas a la izq:
            npos = seqI.index('N')
            seqI = seqI[npos+1:]
            seqI = self.sacNiniciales(nroN, seqI)
        return seqI


    def trimea(self, fasta_in, fasta_out_fh):
        """Toma in FASTA y hace trimming de los extremos en fn de
           el contenido de Ns"""
        for rec in SeqIO.parse(fasta_in, self.__formato_input) :
            laseq = rec.seq.tostring()
            laseq = self.sacNiniciales(10, laseq)
            laseq = self.sacNfinales(10, laseq)
            rec.seq = Seq.Seq(laseq)
            SeqIO.write((rec,), fasta_out_fh, self.__formato_input)
        fasta_out_fh.close()
        return None


    def create_rel(self, XMLin):
        """ Create a dictionary that relate the sequence name
        with the region to mask.

        Returns a dictionary
        """
        bat1 = {}
        b_records = NCBIXML.parse(XMLin)
        for b_record in b_records:
            #print dir(b_record)
            for alin in b_record.alignments:
                #print '='
                #print alin.hit_def
                #print alin.hit_id
                #print alin.title
                for hsp in alin.hsps:
                    qs, qe = hsp.query_start, hsp.query_end
                    if qs > qe:
                        qe, qs = qs, qe
                    bat1.setdefault(b_record.query, set()).add((qs, qe))

        # sort and merge overlapping segments
        for b_record_query in bat1.keys():
            joined_cols = []
            for qs, qe in sorted(list(bat1[b_record_query])):
                if joined_cols:
                    last_qs, last_qe = joined_cols[-1]
                    if last_qe >= qs:
                        joined_cols[-1] = (last_qs, qe)
                        continue
                joined_cols.append((qs, qe))
            bat1[b_record_query] = joined_cols

        return bat1

    def maskseqs(self, ffn, bat1, maskchar='N'):
        """ Take a FASTA file and apply the mask using the
            positions in the dictionary"""
        outseqs_fil = [] # List with SeqRecords
        #outseqs_cut = []
        for record in SeqIO.parse(ffn, self.__formato_input):
            if record.id in bat1:
                # Generate a mutable sequence object to store
                # the sequence with the "mask".
                mutable_seq = record.seq.tomutable()
                #mutable_seq_cut = copy.copy(mutable_seq)
                coords = bat1[record.id]
                for x in coords:
                    mutable_seq[x[0]:x[1]] = maskchar*(x[1]-x[0])
                #mutable_seq_cut = copy.copy(mutable_seq)
                # poner a lo anterior el trimmed.
                # If 15 out of 30 nts in extremes are Ns, mask all 30
                if mutable_seq[:30].count('N')>=15:
                    mutable_seq[:30] = self.__mask*30
                if mutable_seq[-30:].count('N')>=15:
                    mutable_seq[-30:] = self.__mask*30

                seq_rec = SeqRecord(mutable_seq, id=record.id,
                                    description=record.description)
                outseqs_fil.append(seq_rec)
            else:
                # Leave the sequence as found when its name is not
                # in the dictionary.
                outseqs_fil.append(record)
                #outseqs_cut.append(record)
        return outseqs_fil

    def cut(self, fin):
        # TODO: DOC THIS!
        with NamedTemporaryFile(delete=False) as fasta_cut_fh:
            self.trimea(fin, fasta_cut_fh)
            return fasta_cut_fh.name
 

#  blastn -query linker1.txt -subject laurasample1.txt -task blastn-short -evalue 0.0005 -outfmt 5 -out test.xml

    def blastlinkconnector(self, q_seq_fn, s_seq_fn):
        # Make files out of seqs
        with NamedTemporaryFile(delete=False) as blastout_fh:
            blastout_fn = blastout_fh.name
        blastn_cli = blastCL(cmd=self.blastn_path,
                             query=q_seq_fn,
                             subject=s_seq_fn,
                             task="blastn-short",
                             evalue=.00005,
                             outfmt=5, #m 7
                             out=blastout_fn)
        stdout, stderr = blastn_cli()
        return blastout_fn

    def blastregular(self, q_seq_fn, db):
        # Make files out of seqs

        with NamedTemporaryFile(delete=False) as blastout_fh:
            blastout_fn = blastout_fh.name

        blastn_cli = blastCL(cmd=self.blastn_path,
                             query=q_seq_fn,
                             db=db,
                             word_size=11,
                             reward=1,
                             penalty=-5, #-q
                             gapopen=3, #G
                             gapextend=3, #E
                             #dust = "m", #F FIXME
                             #task='blastn',
                             evalue=700, #e o 700?
                             searchsp=1750000000000,
                             outfmt=5, #m 7
                             out=blastout_fn) #o
        stdout, stderr = blastn_cli()

        return blastout_fn


    def apply_filter(self, fin='', db_or_fn='', blast_type='r',
                     mask_char='N', color='vector', bat1_col=''):
        """
        DOC HERE
        """
        #execute = '%s -p blastn -d "%s" -q -5 -G 3 -E 3 -F "m D" -e \
        #25 -Y 1.75e12 -m 7 -i "%s" -o salidablastTMP.xml'%\
        #(self.__blast_path, self.__db_path, fin)
        # TODO: QUE ESTOS PARAMETROS SEAN CFG POR JS!!, POR DEFECTO
        #       SE USE EL DE VECTSCREEN, pero sino sea CFG
        # WHY EFFECTIVESEARCHSP at all?

        if fin == '':
            fin = self.filesfromseqs(self._testseq)

        if blast_type == 'r':
            blastout_fn = self.blastregular(fin, db_or_fn)
        else:
            blastout_fn = self.blastlinkconnector(fin, db_or_fn)

        with open(blastout_fn) as blastout_fh:
            bat1 = self.create_rel(blastout_fh)
        #print open(blastout_fn).read()
        #print bat1
        # make definition dictionary out of bat1

        if not bat1_col:
            bat1_col = {}

        for seqname in bat1:
            bat1_col.setdefault(seqname, {})
            for dnarange in bat1[seqname]:
                bat1_col[seqname][dnarange] = color

        newseqs = self.maskseqs(fin, bat1, mask_char)
        #print newseqs
        # Creates a new temporary file aka filtered
        with NamedTemporaryFile(delete=False) as fasta_out_fh:
            # Write the masked sequence into this temporary file
            SeqIO.write(newseqs, fasta_out_fh, self.__formato_input)
            #print bat1_col
            return fasta_out_fh.name, bat1_col

# si es corrido local, tests

class PaintSeq(object):

    def __init__(self, seq, paint_dict, seq_width=75, isize=''):
        
        # from paint_dict to extended_paint_dict
        self.ex_color_filter = self._get_ext_dict(paint_dict)
        # Add mirna_color
        self.ex_color_filter, self.target_seqs_idxs = \
                              self._getmirna_colors(self.ex_color_filter,
                                             isize)
        self.new_seq_html = self._apply_colors(seq, self.ex_color_filter, seq_width)


    def _get_ext_dict(self, comp_dict):
        """
        Make extended dict out of a compressed dict
        Expected format of comp_dict: {(21, 45): 'connector', (69, 90): 'mirna'}
        """
        ex_color_filter = {}
        for idt in comp_dict:
            for x in range(idt[0]-1, idt[1]):
                ex_color_filter[x] = comp_dict[idt]
        return ex_color_filter

    def _apply_colors(self, seq, ext_paint_dict, seq_width):
        """ """
        new_seq = bytearray()
        for i, x in enumerate(seq):
            if i in ext_paint_dict and (i+1) % seq_width == 0:
                new_seq += '<span class="%s">' % ext_paint_dict[i] + x + '</span><br>'
            elif i in ext_paint_dict and (i+1) % seq_width != 0:
                new_seq += '<span class="%s">' % ext_paint_dict[i] + x + '</span>'
            elif i not in ext_paint_dict and (i+1) % seq_width == 0:
                new_seq += x+'<br>'
            else:
                new_seq += x
        return str(new_seq)

    def _getmirna_colors(self, d, isize):
        """
        Get and colorize target sequences.

        TODO: Make this color independent!!! That is, the user should
        be able to choose color for each kind of filter.
        """

        max_key = max(d)
        connector_color = 0
        vector_color = False
        linker_color = False
        mirna_color = 0
        first_connector_color = False
        for x in range(max_key):
            if x not in d:
                d[x] = 'mirna'
                #print x
                #print 'p3'
            elif not connector_color and not vector_color and d.get(x, '') == 'vector':
                #print 'p1', x
                vector_color = True
            elif not connector_color and d.get(x, '') == 'linker':
                linker_color = True
                # entra en linker_color
                #print 'p2', x, d[x]
            elif not connector_color and linker_color and d.get(x, '') == 'connector' and not first_connector_color:
                vector_color = False
                linker_color = False
                mirna_color = True
                first_connector_color = True
                # genero nuevo array
            elif not vector_color and not linker_color and not connector_color and d.get(x, '') == 'connector' and first_connector_color:
                connector_color += 1
                mirna_color = False
            elif not mirna_color and not vector_color and not linker_color and connector_color and not first_connector_color and d.get(x, '') == 'linker':
                break
            elif connector_color and not vector_color and not linker_color and x not in d:
                d[x] = 'mirna'
                first_connector_color = False
        seqs = []
        if 'linker' in d.values():
            enter_mirna_color = False
            for x in range(max_key):
                # tomar solos los azules
                if enter_mirna_color and d.get(x, '') != 'mirna':
                    print 'NEW: OK', x
                    end_mirna_color = x-1
                    enter_mirna_color = False
                    if end_mirna_color-init_mirna_color>=int(isize):
                        seqs.append((init_mirna_color, end_mirna_color))
                elif d.get(x, '') == 'mirna' and not enter_mirna_color:
                    print 'NEW: OK', x
                    init_mirna_color = x
                    enter_mirna_color = True

        return d, seqs




if __name__ == '__main__':
    #run test
    import doctest
    doctest.testmod()

