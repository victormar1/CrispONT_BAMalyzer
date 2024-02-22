import pysam
from pyfaidx import Fasta

# Chemins
bam_file_path = './barcode19.bam'
reference_fasta_path = './genome.fna'
output_bam_path = './test.bam'
output_unmatched_bam_path = './reads_unmapped.bam'

# Position et chromosome de la condition
positions = [117587768, 117587771]  # Les deux positions où un "C" est requis
chromosome = 'chr7'
nucleotide_target = 'C'

# Ouvrir la séquence de référence
reference = Fasta(reference_fasta_path)

with pysam.AlignmentFile(bam_file_path, "rb") as bam_file, \
     pysam.AlignmentFile(output_bam_path, "wb", template=bam_file) as output_bam, \
     pysam.AlignmentFile(output_unmatched_bam_path, "wb", template=bam_file) as output_unmatched_bam:
    
    for read in bam_file.fetch(chromosome, min(positions)-1, max(positions)):
        reference_positions = read.get_reference_positions(full_length=True)
        sequence = read.query_sequence
        match_conditions = True
        
        for position in positions:
            # Convertir la position en indexation 0-based pour correspondre à Python et à l'indexation de pysam
            position_index = position - 1
            
            # Trouver l'index dans la séquence du read qui correspond à la position du génome
            try:
                read_index = reference_positions.index(position_index)
                if sequence[read_index].upper() != nucleotide_target:
                    match_conditions = False
                    break
            except ValueError:
                # La position n'est pas couverte par le read
                match_conditions = False
                break
        
        if match_conditions:
            output_bam.write(read)
        else:
            output_unmatched_bam.write(read)

reference.close()
