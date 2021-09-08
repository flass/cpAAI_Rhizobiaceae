#!/usr/bin/env python3

import subprocess
import sys
import os
import shutil
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
				label = line.rstrip('\n').split()[0]
				if i==0:
					currlabel = label
				else:
					if label != currlabel:
						raise IndexError("{}\n{}\nlabels (fasta headers) are not orderred the same in input files".format(line, currlabel))
				break # 'for line in fin' loop
			else:
				fout.write(line)
	fout.flush()
	return currlabel

def main(outdir, nflnfmarkgeneseqs=None, nflnfmarkprotseqs=None, nflnfquerygenomes=None, nflnfqueryproteomes=None, nflnfmarkgenealns=None, nflnfmarkprotalns=None, tmpdir='.', nbthreads=1, aligner='mafft', cleanres=False, cleantmp=False, cleanaft=False, reusetmp=False, verbose=False):
	
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
	
	tmpblastdb = path.join(tmpdir, 'blastdb')
	tmpblastout = path.join(tmpdir, 'blastout')
	tmpcombseq = path.join(tmpdir, 'combined_marker_prot_seqs')
	outnuseqdir = path.join(outdir, 'extracted_marker_nucl_seqs')
	outaaseqdir = path.join(outdir, 'extracted_marker_prot_seqs')
	outalndir = path.join(outdir, 'aligned_marker_prot_seqs')
	nfoutconcatprotaln = path.join(outdir, 'concatenated_marker_proteins.aln')
	
	if cleanres:
		print("cleaning: removing previous result files")
		for dirout in [outnuseqdir, outaaseqdir, outalndir]:
			if path.isdir(dirout): shutil.rmtree(dirout)
		if path.isfile(nfoutconcatprotaln): os.remove(nfoutconcatprotaln)
	
	if cleantmp:
		print("cleaning: removing previous temporary files")
		for dirtmp in [tmpblastdb, tmpblastout, tmpcombseq]:
			if path.isdir(dirtmp): shutil.rmtree(dirtmp)
	
	for d in [tmpdir, tmpblastdb, tmpblastout, tmpcombseq, outnuseqdir, outaaseqdir, outalndir]:
		os.makedirs(d, exist_ok=True)
	
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
		nfgenomedb = path.join(tmpblastdb, path.basename(nfquerygenome))
		if not (reusetmp and path.isfile(nfgenomedb)):
#			os.symlink(nfquerygenome, nfgenomedb)
			shutil.copyfile(nfquerygenome, nfgenomedb)
			subprocess.run(['makeblastdb', '-dbtype', 'nucl', '-in', nfgenomedb, '-input_type', 'fasta'])
		
		genomeseqrecs = dict([(r.id, r) for r in SeqIO.parse(nfquerygenome, format='fasta')])
		
		for nfmarkprotseqs in lnfmarkprotseqs:	
			marker = path.splitext(path.basename(nfmarkprotseqs))[0]
			# run blast
			nfblastout = path.join(tmpblastout, genomeid+'_vs_'+marker+'.blastout')
			if not (reusetmp and path.isfile(nfblastout)):
				blastcmd = ['tblastn', '-num_threads', str(nbthreads), '-query', nfmarkprotseqs, '-db', nfgenomedb, '-outfmt', '6', '-out', nfblastout, '-evalue', str(evalthresh)]
				if verbose: print(' '.join(blastcmd))
				subprocess.run(blastcmd)
			
			with open(nfblastout, 'r') as fblastout:
				blastout = [dict(zip(blastoutfmt6stdfields, line.rstrip('\n').split('\t'))) for line in fblastout]
			if len(blastout)>0:
				scores = [float(hit['bitscore']) for hit in blastout]
				maxscore = max(scores)
				besthit = blastout[scores.index(maxscore)]
				hitpos = [besthit[f] for f in ('saccver', 'sstart', 'send')]
				genomehitrec = genomeseqrecs[hitpos[0]]
				if hitpos[1]<hitpos[2]:
					hitbeg = int(hitpos[1]) - 1
					hitend = int(hitpos[2])
				else:
					hitbeg = int(hitpos[2]) - 1
					hitend = int(hitpos[1])
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
		nfoutgeneseq = path.join(outnuseqdir, 'extracted_'+marker+'.fna')
		SeqIO.write(doutgeneseqrec[marker], nfoutgeneseq, format='fasta')
		nfoutprotseq = path.join(outaaseqdir, 'extracted_'+marker+'.faa')
		SeqIO.write(doutprotseqrec[marker], nfoutprotseq, format='fasta')
		
		nfalnout = path.join(outalndir, 'aligned_'+marker+'.faaln')
		if not (reusetmp and path.isfile(nfalnout)):
			if marker in dnfmarkprotalns:
				print("aligning extracted protein sequences for marker {} to input reference alignment".format(marker))
				if aligner=='clustalo':
					aligncmd = ['clustalo', '--in', nfoutprotseq, '--profile1', dnfmarkprotalns[marker], '-t', 'Protein', '--threads', str(nbthreads), '--out', nfalnout, '--infmt', 'fasta', '--outfmt', 'fasta']
				elif aligner=='mafft':
					aligncmd = ['mafft', '--quiet', '--amino', '--inputorder', '--thread', str(nbthreads), '--add', nfoutprotseq, dnfmarkprotalns[marker]]
			else:
				print("aligning extracted protein sequences for marker {} together with input reference sequences".format(marker))
				nftmpcombseq = path.join(tmpcombseq, marker+'.faa')
				shutil.copyfile(nfmarkprotseqs, nftmpcombseq)
				with open(nftmpcombseq, 'a') as ftmpcombseq:
					with open(nfoutprotseq, 'r') as foutprotseq:
						ftmpcombseq.write(foutprotseq.read())
				if aligner=='clustalo':
					aligncmd = ['clustalo', '--in', nftmpcombseq, '-t', 'Protein', '--threads', str(nbthreads), '--out', nfalnout, '--infmt', 'fasta', '--outfmt', 'fasta']
				elif aligner=='mafft':
					aligncmd = ['mafft', '--quiet', '--amino', '--inputorder', '--thread', str(nbthreads), nftmpcombseq]
			if verbose: print(' '.join(aligncmd))
			if aligner=='clustalo':
				subprocess.run(aligncmd)
			elif aligner=='mafft':
				with open(nfalnout, 'w') as falnout:
					subprocess.run(aligncmd, stdout=falnout)
	
	# concatenate protein alignments
	lfinhandles = [open(path.join(outalndir, bnfalnout), 'r') for bnfalnout in os.listdir(outalndir)]
	with open(nfoutconcatprotaln, 'w') as foutconcatprotaln:
		print("concatenating the marker protein alignments")
		currlabel = None
		nextlabel = iterOneLabel(lfinhandles, foutconcatprotaln, currlabel)
		n = 1
		while nextlabel != currlabel:
			currlabel = nextlabel
			nextlabel = iterOneLabel(lfinhandles, foutconcatprotaln, currlabel)
			if verbose: sys.stdout.write("\r{} {}".format(n, nextlabel.rstrip('\n')))
			n += 1

		for fin in lfinhandles:
			fin.close()
		if verbose: sys.stdout.write("\n")

	if cleanaft:
		print("cleaning: removing temporary files from this run")
		for dirtmp in [tmpblastdb, tmpblastout, tmpcombseq]:
			os.rmdir(dirtmp)


def usage():
	s = "genome2cpAAI.py {-g/--list_marker_gene_seqs listfilepaths|-p/--list_marker_prot_seqs listfilepaths} {-q/--list_query_genomes listfilepaths|-Q/--list_query_proteomes listfilepaths} [{-a/--list_marker_gene_alns listfilepaths|-A/--list_marker_prot_alns listfilepaths}] [--threads int(default: 1)] [--tmp_dir folderpath(default: current directory)] [--cleanres] [--cleantmp] [--cleanaft] [--cleanall] [--aligner {mafft|clustalo}(default: mafft)']"
	return s

if __name__=="__main__":
	
	opts, args = getopt.getopt(sys.argv[1:], 'p:g:a:A:q:Q:o:t:hv', 
							   ['list_marker_gene_seqs=', 'list_marker_prot_seqs=', \
								'list_marker_gene_alns=', 'list_marker_prot_alns=', \
								'list_query_genomes=', 'list_query_proteomes=', \
								'ouput_dir=', 'aligner=', 'tmp_dir=', 'reuse_prevtmp', \
								'clean_prevres', 'clean_prevtmp', 'clean_after', 'clean_all', \
								'threads=', 'verbose', 'help'])
	dopts = dict(opts)
	if ('-h' in dopts) or ('--help' in dopts):
		print(usage())
		sys.exit(0)
	
	nflnfmarkgeneseqs = dopts.get('-g', dopts.get('--list_marker_gene_seqs'))
	if nflnfmarkgeneseqs:
		raise ValueError("argument -g/--list_marker_gene_seqs are not currently supported")
	nflnfmarkprotseqs = dopts.get('-p', dopts.get('--list_marker_prot_seqs'))
	if not (nflnfmarkgeneseqs or nflnfmarkprotseqs):
		raise ValueError("missing mandatory argument; must provide a value to either arguments: {-g/--list_marker_gene_seqs|-p/--list_marker_prot_seqs}")

	nflnfmarkgenealns = dopts.get('-a', dopts.get('--list_marker_gene_alns'))
	if nflnfmarkgenealns:
		raise ValueError("argument -a/--list_marker_gene_alns are not currently supported")
	nflnfmarkprotalns = dopts.get('-A', dopts.get('--list_marker_prot_alns'))
	
	nflnfquerygenomes = dopts.get('-q', dopts.get('--list_query_genomes'))
	nflnfqueryproteomes = dopts.get('-Q', dopts.get('--list_query_proteomes'))
	if nflnfqueryproteomes:
		raise ValueError("argument -Q/--list_query_proteomes are not currently supported")
	if not (nflnfquerygenomes or nflnfqueryproteomes):
		raise ValueError("missing mandatory argument; must provide a value to either arguments: {-q/--list_query_genomes|-Q/--list_query_proteomes}")
	
	nbthreads = int(dopts.get('-t', dopts.get('--threads', 1)))
	
	outdir = dopts.get('-o', dopts.get('--ouput_dir'))
	if not outdir:
		raise ValueError("missing mandatory argument; must provide a value to argument -o/--output_dir")
	
	tmpdir = dopts.get('--tmp_dir', '.')
	
	aligner = dopts.get('--aligner', 'mafft')
	
	cleanres = ('--clean_prevres' in dopts)
	reusetmp = ('--reuse_prevtmp' in dopts)
	cleantmp = ('--clean_prevtmp' in dopts)
	if reusetmp and cleantmp: raise ValueError("cannot combine options --reuse_prevtmp and --clean_prevtmp")
	cleanaft = ('--clean_after' in dopts)
	if ('--clean_all' in dopts):
		cleanres = cleantmp = cleanaft = True
	
	verbose = ('--verbose' in dopts) or ('-v' in dopts) 

	main(outdir, nflnfmarkgeneseqs=nflnfmarkgeneseqs, nflnfmarkprotseqs=nflnfmarkprotseqs, nflnfquerygenomes=nflnfquerygenomes, nflnfqueryproteomes=nflnfqueryproteomes, nflnfmarkgenealns=nflnfmarkgenealns, nflnfmarkprotalns=nflnfmarkprotalns, tmpdir=tmpdir, nbthreads=nbthreads, aligner=aligner, cleanres=cleanres, cleantmp=cleantmp, cleanaft=cleanaft, reusetmp=reusetmp, verbose=verbose)
