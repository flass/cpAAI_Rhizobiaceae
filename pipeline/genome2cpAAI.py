#!/usr/bin python3

import subprocess
import sys
import os
from os import path
import getopt
from Bio import SeqIO

blastoutfmt6stdfields = ['qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
evalthresh = 0.00001
transtable = 11

# function taken from concatAlnSameLabels
def iterOneLabel(lfinhandles, fout, currlabel):
	if not (currlabel is None):
		fout.write(currlabel+'\n')
	for i, fin in enumerate(lfinhandles):
		for line in fin:
			if line.startswith('>'):
				if i==0:
					currlabel = line.rstrip('\n').split()[0]
				else:
					if line.rstrip('\n').split()[0] != currlabel:
						print line
						print currlabel
						raise IndexError, "labels (fasta headers) are not orderred the same in input files"
				break # 'for line in fin' loop
			else:
				fout.write(line)
	fout.flush()
	return currlabel

def main(outdir, nflnfmarkgeneseqs=None, nflnfmarkprotseqs=None, nflnfquerygenomes=None, nflnfqueryproteomes=None, nflnfmarkgenealns=None, nflnfmarkprotalns=None, tmpdir='.', nbthreads=1, cleantmp=False):
	
	if nflnfquerygenomes:
		with open(nflnfquerygenomes, 'r') as flnfquerygenomes:
			lnfquerygenomes = [line.rstrip('\n') for line in flnfquerygenomes]
	elif nflnfqueryproteomes:
		raise Warning("not supported yet")
		with open(nflnfqueryproteomes, 'r') as flnfqueryproteomes:
			lnfqueryproteomes = [line.rstrip('\n') for line in flnfqueryproteomes]
	else:
		raise ValueError("need to provide input list of either genomes or proteomes")
		
	if nflnfmarkgeneseqs:
		raise Warning("not supported yet")
		with open(nflnfmarkgeneseqs, 'r') as flnfmarkgeneseqs:
			lnfmarkgeneseqs = [line.rstrip('\n') for line in flnfmarkgeneseqs]
	elif nflnfmarkprotseqs:
		with open(nflnfmarkprotseqs, 'r') as flnfmarkprotseqs:
			lnfmarkprotseqs = [line.rstrip('\n') for line in flnfmarkprotseqs]
	else:
		raise ValueError("need to provide input list of either marker gene or marker protein sequence Fasta files (one file per marker locus/protein, multiple sequences in each)")
	
	dnfmarkgenealns = {}
	dnfmarkprotalns = {}
	if nflnfmarkgenealns:
		raise Warning("not supported yet")
		if not nflnfmarkgeneseqs:
			raise ValueError("need to provide input gene marker sequences consistent with gene marker alignments")
		with open(nflnfmarkgenealns, 'r') as flnfmarkgenealns:
			lnfmarkgenealns = [line.rstrip('\n') for line in flnfmarkgenealns]
			dnfmarkgenealns = dict([(path.splitext(path.basename(nfmarkgenealns))[0], nfmarkgenealns) for nfmarkgenealns in lnfmarkgenealns])
	elif nflnfmarkprotalns:
		if not nflnfmarkprotseqs:
			raise ValueError("need to provide input protein marker sequences consistent with protein marker alignments")
		with open(nflnfmarkprotalns, 'r') as flnfmarkprotalns:
			lnfmarkprotalns = [line.rstrip('\n') for line in flnfmarkprotalns]
			dnfmarkprotalns = dict([(path.splitext(path.basename(nfmarkprotaln))[0], nfmarkprotaln) for nfmarkprotaln in lnfmarkprotalns])
	else:
		print("no marker gene/protein alignment provided, will have to align marker gene/protein sequences from scratch together with extracted input")
	
	os.makedirs(tmpdir, exist_ok=True)
	tmpblastdb = path.join(tmpdir, 'blastdb')
	os.makedirs(tmpblastdb, exist_ok=True)
	tmpblastout = path.join(tmpdir, 'blastout')
	os.makedirs(tmpblastout, exist_ok=True)
	tmpcombseq = path.join(tmpdir, 'combined_marker_prot_seqs')
	os.makedirs(tmpcombseq, exist_ok=True)
	outseqdir = path.join(outdir, 'extracted_marker_prot_seqs')
	os.makedirs(outseqdir, exist_ok=True)
	outalndir = path.join(outdir, 'aligned_marker_prot_seqs')
	os.makedirs(outalndir, exist_ok=True)
	
	# dict will store extracted sequence records of marker genes
	doutgeneseqrec = {}
	doutprotseqrec = {}
	for nfmarkprotseqs in lnfmarkprotseqs:
		marker = path.splitext(path.basename(nfmarkprotseqs))[0]
		doutgeneseqrec[marker] = []
		doutprotseqrec[marker] = []
	
	for nfquerygenome in lnfquerygenomes:
#		qgenomeseqrecs = SeqIO.parse(nfquerygenomes, format='fasta')
		genomeid = path.splitext(path.basename(nfquerygenome))[0]
		# make blast db
	    markprotseqdb = path.join([tmpblastdb, path.basename(nfmarkprotseqs)])
	    os.symlink(nfmarkprotseqs, markprotseqdb)
		subprocess.run(['makeblastdb', '-dbtype', 'prot', '-in', markprotseqdb, '-input_type', 'fasta'])
		genomedb = path.join([tmpblastdb, path.basename(nfquerygenome)])
	    os.symlink(nfquerygenome, genomedb)
		subprocess.run(['makeblastdb', '-dbtype', 'nucl', '-in', genomedb, '-input_type', 'fasta'])
		
		genomeseqrecs = dict([(r.id, r) for r in SeqIO.parse(nfquerygenome, format='fasta')])
		
		for nfmarkprotseqs in lnfmarkprotseqs:	
			marker = path.splitext(path.basename(nfmarkprotseqs))[0]
			# run blast
			nfblastout = path.join([tmpblastout, genomeid+'_vs_'+marker+'.blastout'])
			subprocess.run(['tblastn', '-num_threads', str(nbthreads), '-query', nfquerygenome, '-db', nfmarkprotseqs, '-outfmt', '6', '-out', nfblastout, '-evalue', str(evalthresh)])
			
			with open(nfblastout, 'r') as fblastout:
				blastout = [dict(zip(blastoutfmt6stdfields, line.rstrip('\n').split('\t'))) for line in flastout]
			if len(blastout)>0:
				scores = [hit['bitscore'] for hit in blastout]
				maxscore = max(scores)
				besthit = blastout[scores.index(maxscore)]
				hitpos = [besthit[f] for f in ('saccver', 'sbeg', 'send', 'sframe')]
				genomehitrec = genomeseqrecs[hitpos[0]]
				if hitpos[1]<hitpos[2]:
					hitbeg = hitpos[1] - 1
					hitend = hitpos[2]
				else:
					hitbeg = hitpos[2] - 1
					hitend = hitpos[1]
				hitseqrec = genomehitrec[hitbeg:hitend]
#				hitseqrec.id = genomeid+'_'+marker
				hitseqrec.id = genomeid
				hitseqrec.name += '_'+marker
				hitseqrec.description += ' [{}..{}]'.format(hitbeg,hitend)
				if hitpos[1]>hitpos[2]:
					hitseqrec.seq = hitseqrec.seq.reverse_complement()
					hitseqrec.description += ' (reverse complement)'
				doutgeneseqrec[marker].append(hitseqrec)
				# translate to protein
				protseqrec = hitseqrec.translate(table=transtable)
				protseqrec.id = hitseqrec.id
				protseqrec.name = hitseqrec.name
				protseqrec.description = 'translation of '+hitseqrec.description
				doutprotseqrec[marker].append(protseqrec)


	# for each marker write combined extracted sequences from all genomes to output files
	dfoutseq = {}
	for nfmarkprotseqs in lnfmarkprotseqs:
		marker = path.splitext(path.basename(nfmarkprotseqs))[0]
		nfoutgeneseq = path.join([outseqdir, 'extracted_'+marker+'.fna'])
		SeqIO.write(doutgeneseqrec[marker], nfoutgeneseq, format='fasta')
		nfoutprotseq = path.join([outseqdir, 'extracted_'+marker+'.faa'])
		SeqIO.write(doutprotseqrec[marker], nfoutprotseq, format='fasta')
		
		nfalnout = path.join([outalndir, 'aligned_'+marker+'.faaln'])
		if marker in dnfmarkprotalns:
			print("aligning extracted protein sequences for marker {} to input reference alignment".format(marker))
			subprocess.run(['clustalo', '-i', doutprotseqrec[marker], '--profile1', dnfmarkprotalns[marker], '-t', 'Protein', '--threads', str(nbthreads), '-o', nfalnout, '-infmt', 'fasta', '-outfmt', 'fasta'])
		else:
			print("aligning extracted protein sequences for marker {} to input reference sequences".format(marker))
			nftmpcombseq = path.join([tmpcombseq, marker+'.faa'])
			os.copy(nfmarkprotseqs, nftmpcombseq)
			with open(nftmpcombseq, 'a') as ftmpcombseq:
				with open(nfoutprotseqm, 'r') as foutprotseq:
					ftmpcombseq.write(foutprotseq.read())
			subprocess.run(['clustalo', '-i', ftmpcombseq, '-t', 'Protein', '--threads', str(nbthreads), '-o', nfalnout, '-infmt', 'fasta', '-outfmt', 'fasta'])
	
	# concatenate protein alignments
	lfinhandles = [open(path.join(outalndir, bnfalnout), 'r') for bnfalnout in os.listdir(outalndir)]
	with open(path.join(outdir, 'concatenated_marker_proteins.aln'), 'w') as foutconcatprotaln:
		currlabel = None
		nextlabel = iterOneLabel(lfinhandles, foutconcatprotaln, currlabel)
		n = 1
		while nextlabel != currlabel:
			currlabel = nextlabel
			nextlabel = iterOneLabel(lfinhandles, foutconcatprotaln, currlabel)
			sys.stdout.write("\r{} {%s}".format(n, nextlabel.rstrip('\n')))
			n += 1

		for fin in lfinhandles:
			fin.close()
	
	if cleantmp:
		print("cleaning: removing temporary files")
		for dirtmp in [tmpblastdb, tmpblastout, tmpcombseq]:
			os.rmdir(dirtmp)


def usage():
	s = "{} {-g/--list_marker_gene_seqs listfilepaths|-p/--list_marker_prot_seqs listfilepaths} {-q/--list_query_genomes listfilepaths|-Q/--list_query_proteomes listfilepaths} [{-a/--list_marker_gene_alns listfilepaths|-A/--list_marker_prot_alns listfilepaths}] [--threads int(default: 1)] [--tmp_dir folderpath(default: current directory)] [--cleantmp]".format(sys.argv[0])
	return s

if __name__=="__main__":
	
	opts, args = getopt.getopt(sys.argv[1:], 'p:g:a:A:q:Q:t:ho:', 
							   ['list_marker_gene_seqs=', 'list_marker_prot_seqs=', \
								'list_marker_gene_alns=', 'list_marker_prot_alns=', \
								'list_query_genomes=', 'list_query_proteomes=', \
								'ouput_dir=', 'tmp_dir=', 'cleantmp', \
								'threads=', 'help'])
	
	if ('-h' in opts) or ('--help' in opts):
		print(usage())
		sys.exit(0)
	
	nflnfmarkgeneseqs = opts.get('-g', opts.get('--list_marker_gene_seqs'))
	if nflnfmarkgenealns:
		raise ValueError("argument -g/--list_marker_gene_seqs are not currently supported")
	nflnfmarkprotseqs = opts.get('-p', opts.get('--list_marker_prot_seqs'))
	if not (nflnfmarkgeneseqs or nflnfmarkprotseqs):
		raise ValueError("missing mandatory argument; must provide either arguments: {-g/--list_marker_gene_seqs|-p/--list_marker_prot_seqs}")

	nflnfmarkgenealns = opts.get('-a', opts.get('--list_marker_gene_alns'))
	if nflnfmarkgenealns:
		raise ValueError("argument -a/--list_marker_gene_alns are not currently supported")
	nflnfmarkprotalns = opts.get('-A', opts.get('--list_marker_prot_alns'))
	
	nflnfquerygenomes = opts.get('-q', opts.get('--list_query_genomes'))
	if nflnfquerygenomes:
		raise ValueError("argument -q/--list_query_genomes are not currently supported")
	nflnfqueryproteomes = opts.get('-Q', opts.get('--list_query_proteomes'))
	if not (nflnfquerygenomes or nflnfqueryproteomes):
		raise ValueError("missing mandatory argument; must provide either arguments: {-q/--list_query_genomes|-Q/--list_query_proteomes}")
	
	nbthreads = int(opts.get('-t', opts.get('--threads', 1)))
	
	outdir = opts.get('-o', opts.get('--ouput_dir'))
	tmpdir = opts.get('--tmp_dir', '.')
	cleantmp = opts.get('--cleantmp', False)

	main(outdir, nflnfmarkgeneseqs=nflnfmarkgeneseqs, nflnfmarkprotseqs=nflnfmarkprotseqs, nflnfquerygenomes=nflnfquerygenomes, nflnfqueryproteomes=nflnfqueryproteomes, nflnfmarkgenealns=nflnfmarkgenealns, nflnfmarkprotalns=nflnfmarkprotalns, tmpdir=tmpdir, nbthreads=nbthreads, cleantmp=cleantmp)
