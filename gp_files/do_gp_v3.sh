#!/bin/bash

pepsi_path=/home/malab/Desktop/Pepsi-SAXS-Linux/Pepsi-SAXS
structures=$1
exp_path=$2
theta=$3  # Currently unused, unless iBME is activated
gl=$4 #Optional

gpdoc=$5

grep -v "#" $2 | awk '{print $1}' > qvals.txt


# Check if structure directory exists
if [[ ! -d "$structures" ]]; then
    echo "Error: structures not found: $structures"
    exit 1
fi

# Check if experimental data file exists
if [[ ! -f "$exp_path" ]]; then
    echo "Error: experimental data not found: $exp_path"
    exit 1
fi

# Get sorted list of structure files named mm016_*.pdb
structure_files=( $(find "$structures" -name "mm016_*.pdb" | sort) )
ens_size=${#structure_files[@]}

# Choose mode: single grid point or all
if [[ -n "$gl" ]]; then
    # ----------- SINGLE GRID POINT MODE -----------
    total_lines=$(grep -v "#" $5 | wc -l)
    if (( gl < 1 || gl > total_lines )); then
        echo "Error: Grid line $gl is out of range (1 to $total_lines)"
        exit 1
    fi
    grid_lines=($gl)
else
    # ----------- ALL GRID POINTS MODE -------------
    grid_lines=( $(seq 1 $(grep -v "#" $5 | wc -l)) )
fi

# Loop over selected grid point lines
for gl in "${grid_lines[@]}"; do
    grid_point=$(grep -v "#" $5 | sed -n "${gl}p" | awk '{print $1}')
    dro=$(grep -v "#" $5 | sed -n "${gl}p" | awk '{print $2}')
    r0=$(grep -v "#" $5 | sed -n "${gl}p" | awk '{print $3}')

    echo "Running Pepsi-SAXS for grid point $grid_point (line $gl)..."
    mkdir -p GP$grid_point
  
    > GP$grid_point/calc_saxs.txt   # clear old data
    # Loop over structure files
    for i in "${!structure_files[@]}"; do
	echo "Processing structure $i: ${structure_files[$i]}"
        pdb_file="${structure_files[$i]}"
        $pepsi_path "$pdb_file" "$exp_path" -o GP$grid_point/saxs$i.dat \
            -cst --cstFactor 0 --I0 1.0 --dro $dro \
            --r0_min_factor $r0 --r0_max_factor $r0 --r0_N 1
       # Extract SAXS intensity column (q, I(q), etc.) — store only I(q)
        intensities=$(awk '!/^#/ {printf "%s ", $4}' GP$grid_point/saxs$i.dat)
	echo "$i $intensities" >> GP$grid_point/calc_saxs.txt
            
            

        rm GP$grid_point/saxs$i.dat

        # Extract Rg only for grid_point == 0
        if [ "$grid_point" -eq 0 ]; then
            grep "Radius of gyration of the envelope" GP$grid_point/saxs$i.log | \
                awk '{print $8}' >> GP$grid_point/Rg_env.dat
        fi

        rm GP$grid_point/saxs$i.log
    done > GP$grid_point/logPEPSI

    # Move any generated iBME files (if run separately)
    #mv gp${grid_point}_* GP$grid_point/ 2>/dev/null
    #mv gp${grid_point}.log GP$grid_point/ 2>/dev/null
done
