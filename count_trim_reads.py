import re
import gzip
import glob

# https://github.com/lh3/readfq/blob/master/readfq.py
def readfq(fp): 
    # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

files = glob.glob('/scratch/drabosky_flux/sosi/AWT_delimit/trim_reads/*R1*')

for file1 in files:
	file2 = re.sub('R1', 'R2', file1)
	fileun = re.sub('R1', 'unpaired', file1)

	totlen = 0
	for f in [file1, file2, fileun]:
		# print(f)
		f = gzip.open(f, 'r')
		for idname, seq, qual in readfq(f):
			totlen += len(seq)
		f.close()

	sample = re.sub('^.*\/', '', file1)
	sample = re.sub('_R1.*', '', sample)

	print('%s,%s' % (sample, totlen))
