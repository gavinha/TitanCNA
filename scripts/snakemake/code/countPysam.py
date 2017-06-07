# requires python3
# requires pysam-0.11.2.1
import sys
import pysam

chromToUse = sys.argv[1]  # 0 for all chromosomes
norm_hetpsns = sys.argv[2]
bam_file = sys.argv[3]
#ref_file = sys.argv[4]
base_quality = int(sys.argv[4])
map_quality = int(sys.argv[5])
vcf_quality = int(sys.argv[6])
positions = {}

#  add (position,depth) from the normal hetpositions input file to a dictionary of lists
#    indexed by chromosome
for line in open(norm_hetpsns):
    if not line.strip().startswith("#"):
        chrom = line.split()[0] 
        if chrom == chromToUse or chromToUse == 0:
            position = int(line.strip().split()[1])
            ref_base = line.strip().split()[3]
            nref_base = line.strip().split()[4]
            qual = line.strip().split()[5]
            depth = line.split()[7].split(';')[0].replace('DP=', '') 
            position_data = position, depth, ref_base, nref_base, qual
            if chrom not in positions:
                positions[chrom] = []
            positions[chrom].append(position_data)

sample = pysam.AlignmentFile(bam_file)
#reference = pysam.FastaFile(ref_file)
## print header ##
print ("Chr\tPosition\tRef\tRefCount\tNref\tNrefCount\tNormQuality")

for chrom in positions:
  i = 0
  for position_data in positions[chrom]:
  	position = int(position_data[0])
  	result = str(chrom) + "\t" + str(position)
  	ref_base = position_data[2]
  	nref_base = position_data[3]
  	qual = float(position_data[4])
  	if qual >= vcf_quality and qual != None:
  		_p = sample.pileup(reference=chrom, start=position, end=position + 1)
  		bases = list()
  		for p in _p:
  			if p.reference_pos == position:
  				for r in p.pileups:
  					if not r.is_del and not r.is_refskip:
  						base = r.alignment.query_sequence[r.query_position-1]
  						mapq = r.alignment.mapping_quality
  						baseq = r.alignment.query_qualities[r.query_position-1]
  						if mapq >= map_quality and baseq >= base_quality:
  							bases.append(base)
  		ref_count = 0
  		depth = 0
  		for base in bases:
  			depth += 1
  			if base == ref_base:
  				ref_count += 1
  		alt_count = depth - ref_count
  		
  		result += "\t" + ref_base + "\t" + str(ref_count) + "\t" + nref_base + "\t" + str(alt_count) + '\t' + str(qual)
  		print(result)
  		i += 1
						
        
