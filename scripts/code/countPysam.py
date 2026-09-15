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
print("Chr\tPosition\tRef\tRefCount\tNref\tNrefCount\tNormQuality")

for chrom in positions:
    # Build SNP lookup dictionary, filtering by vcf_quality
    snp_lookup = {}
    for position_data in positions[chrom]:
        pos_1based = int(position_data[0])
        ref_base   = position_data[2]
        nref_base  = position_data[3]
        qual       = float(position_data[4])
        if qual >= vcf_quality:
            snp_lookup[pos_1based] = (ref_base, nref_base, qual)

    # Skip chromosome if no positions pass quality filter
    if not snp_lookup:
        continue

    # Create set for O(1) membership checks and find bounds
    snp_positions = set(snp_lookup.keys())
    min_pos = min(snp_positions)
    max_pos = max(snp_positions)

    # Collect pileup data for all SNP positions with coverage
    # Key: position (1-based), Value: list of bases passing quality filters
    pileup_data = {}

    # Single pileup call per chromosome spanning all SNPs
    for pileupcolumn in sample.pileup(
        reference=chrom,
        start=min_pos - 1,  # Convert 1-based to 0-based (pysam is 0-based)
        end=max_pos,        # Half-open interval: covers up to max_pos (1-based)
        truncate=False,
        stepper="samtools"
    ):
        # Convert pysam 0-based position back to VCF 1-based
        vcf_pos = pileupcolumn.pos + 1

        # Skip positions that are not in our SNP set
        if vcf_pos not in snp_positions:
            continue

        # Process reads at this SNP position
        bases = []
        for r in pileupcolumn.pileups:
            if not r.is_del and not r.is_refskip:
                base  = r.alignment.query_sequence[r.query_position]
                mapq  = r.alignment.mapping_quality
                baseq = r.alignment.query_qualities[r.query_position]
                if mapq >= map_quality and baseq >= base_quality:
                    bases.append(base)

        pileup_data[vcf_pos] = bases

    # Iterate through SNP positions in original order and print output
    # Positions failing vcf_quality are skipped; zero-coverage positions are still written
    for position_data in positions[chrom]:
        position  = int(position_data[0])
        ref_base  = position_data[2]
        nref_base = position_data[3]
        qual      = float(position_data[4])

        if qual >= vcf_quality:
            bases     = pileup_data.get(position, [])
            ref_count = bases.count(ref_base)
            alt_count = len(bases) - ref_count

            result = (f"{chrom}\t{position}\t{ref_base}\t{ref_count}"
                      f"\t{nref_base}\t{alt_count}\t{qual}")
            print(result)