import pandas as pd
import subprocess
import os

# Species list
species = ['bonobo', 'human', 'gorilla', 'macaque', 'marmoset', 'chimpanzee']

# Minimum fragments threshold (adjust as needed)
min_frags_threshold = 1000


# Create all commands first
commands = []
for _, row in subsample_summary.iterrows():
    celltype = row['celltype']
    celltype = celltype.lower()
    ncells = int(row['ncells'])
    nfrags = int(row['nfragments'])
    frags_per_cell = int(row['avg_fragments'])
    min_frags_threshold = frags_per_cell*0.75
    
    output_dir = f"subsampled_fragments/minimum_{celltype}_{nfrags}"
    os.makedirs(output_dir, exist_ok=True)
    
    for sp in species:
        input_file = f"fragment_files/{celltype}/{sp}_{celltype}.fragments.tsv.gz"
        output_file = f"{output_dir}/{sp}_{celltype}.fragments_sub.tsv.gz"
        
        if os.path.exists(input_file):
            cmd = f"bash subsample_fragments.sh -i {input_file} -o {output_file} -c {ncells} -p {frags_per_cell} -m {min_frags_threshold}"
            cmd = f"bash subsample_fragments.sh -i {input_file} -o {output_file} -f {nfrags} -m {min_frags_threshold}"
            commands.append(cmd)

# Write commands to file
with open('subsample_commands.txt', 'w') as f:
    f.write('\n'.join(commands))

print(f"📝 Generated {len(commands)} commands")
print("🚀 Running in parallel with GNU parallel...")

# Run in parallel (adjust -j for number of parallel jobs)
subprocess.run([
    'parallel',
    '-j', '8',  # 8 parallel jobs
    '--eta',    # Show progress
    '--bar',    # Progress bar
    ':::', 
] + commands, check=True)

print("\n🎉 All subsampling complete!")