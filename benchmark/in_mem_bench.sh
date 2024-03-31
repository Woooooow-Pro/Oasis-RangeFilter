#!/bin/bash

# Reference: https://github.com/Eric-R-Knorr/Proteus/blob/master/bench/in-memory.sh

# Exit if a command fails
set -o errexit

# Throw error when accessing an unset variable
set -o nounset

# Enable debug mode with TRACE=1 ./bench.sh
if [[ "${TRACE-0}" == "1" ]]; then
    set -o xtrace
fi

# Change to directory of this script
cd "$(dirname "$0")"

###############################################
###########   WORKLOAD PARAMETERS   ###########
###############################################

## 1. Key Length (in bits)
#######################################
# For integer workloads, this value MUST be 64

## 2. Number of Keys
## 3. Number of Queries
#######################################

## [LINKED]
## 4. Key Distribution
## 5. Query Distribution
#######################################
# For SOSD / Domain workloads, please download the relevant datasets and follow the instructions in the README
# SOSD and Domain Key and Query Distributions must be specified together, e.g. ksosd_books + qsosd_books
# Integer Key Distributions: kuniform, knormal, ksosd_books, ksosd_fb
# Integer Query Distributions: quniform, qnormal, qsosd_books, qsosd_fb

## [LINKED]
## 6. Minimum Range Size
## 7. Maximum Range Size
## 8. Point / Range Query Ratio
#######################################
# Min Range Size must be >= 2
# Point / Range Query Ratio must be >= 0.0 && <= 1.0
# Use 0.0 for all range queries and 1.0 for all point queries
# remaining point queries are assigned to the uniform distribution.

## 9. Positive / Negative Query Ratio
#######################################
# Must be >= 0.0 && <= 1.0
# Use 0.0 for all negative queries and 1.0 for all positive queries

keylen_arr=(64)
nkeys_arr=(10000000)
nqrys_arr=(1000000)
kdist_arr=("ksosd_fb" "ksosd_books" "kuniform" "knormal")
qdist_arr=("qsosd_fb" "qsosd_books" "quniform" "quniform")
minrange_arr=(2 2 2 2)
maxrange_arr=(2 32 512 1024)
pqratio_arr=(1.0 0.0 0.0 0.5)
pnratio_arr=(0.0)

###############################################
############   FILTER PARAMETERS   ############
###############################################

## Filter Bits-per-Key
#######################################
# Must be a positive real number
membudg_arr=(8 10 12 14 16 18)

## Block Size
#######################################
block_sizes=(150)

##############################################################################################

REPO_DIR="$(pwd)"
SOSD_DIR="$REPO_DIR/workloads/SOSD/"
WORKL_BIN="$REPO_DIR/../build/bin/workload_gen"
EXP_BIN="$REPO_DIR/../build/bin//in_mem_bench"

# Create result folder
res_folder="$REPO_DIR/in_mem_result/$1"
if [ ! -e "$res_folder" ]; then
    mkdir -p "$res_folder"
else
    echo "Result folder already exists! Please try again with a valid result folder name."
    exit -1
fi

SCRATCH_DIR="$REPO_DIR/in_mem_result/tmp"
DATA_DIR="$REPO_DIR/in_mem_data"
mkdir -p "$DATA_DIR"
mkdir -p "$SCRATCH_DIR"

# Create file that indexes the information of the generated experiment result files
index_file="$res_folder/index.txt"
touch $index_file

# Create result csv
res_csv="$res_folder/results.csv"
touch $res_csv

fetch_data() {
    path="$DATA_DIR/${nkeys// /_}/${nqrys// /_}/${minrange// /_}/${maxrange// /_}/${kdist// /_}/${qdist// /_}/${pqratio// /_}/${pnratio// /_}"

    if [ ! -e "$path/my_data" ]; then
        echo "Generating data in $path"
        mkdir -p "$path" && cd "$path"
        $WORKL_BIN "$SOSD_DIR" "$nkeys" "$nqrys" "$minrange" "$maxrange" "$kdist" "$qdist" "$pqratio" "$pnratio"
    else
        echo "Copying data from $path"
    fi

    cp -r "$path/my_data" $expdir/my_data
}

experiment() {
    cd "$expdir"
    touch ./experiment_result
    echo -e "### FILE INFO ###" >>./experiment_result
    echo -e "\tResult Folder:\t$res_folder" >>./experiment_result
    echo -e "\tFile Counter:\t$filecnt\n\n" >>./experiment_result

    echo -e "### BEGIN STANDALONE EXPERIMENT ###" >>./experiment_result

    echo -e "\n\n### BEGIN EXPERIMENT DESCRIPTION ###" >>./experiment_result
    echo -e "\tFilter Name:\t$filter" >>./experiment_result
    echo -e "\tNumber of Keys:\t$nkeys" >>./experiment_result
    echo -e "\tKey Length:\t$keylen" >>./experiment_result
    echo -e "\tNumber of Queries:\t$nqrys" >>./experiment_result
    echo -e "\tKey Distribution:\t$kdist" >>./experiment_result
    echo -e "\tQuery Distribution:\t$qdist" >>./experiment_result
    echo -e "\tMinimum Range Size:\t$minrange" >>./experiment_result
    echo -e "\tMaximum Range Size:\t$maxrange" >>./experiment_result
    echo -e "\tPoint / Range Query Ratio:\t$pqratio" >>./experiment_result
    echo -e "\tPositive / Negative Query ratio:\t$pnratio" >>./experiment_result

    fetch_data
    cd "$expdir"

    # Print experiment parameters to index file and result CSV
    out_file="$filter-$filecnt"
    echo -e "############# $out_file #############" >>$index_file
    echo -e "NKeys: $nkeys; NQueries: $nqrys; KLen: $keylen" >>$index_file
    echo -e "KDist: $kdist; QDist: $qdist" >>$index_file
    echo -e "MinRange: $minrange; MaxRange: $maxrange" >>$index_file
    echo -e "P-Q ratio: $pqratio; +/- Q ratio: $pnratio" >>$index_file

    config="$filter-${nkeys// /_}-${keylen// /_}-${nqrys// /_}-${minrange// /_}-${maxrange// /_}-${pqratio// /_}-${kdist// /_}-${qdist// /_}-${pnratio// /_}"
    if [ $filter = "OasisPlus" ]; then
        config+="-${membudg}-${block_size}"
        echo -e "BPK: $membudg; BlockSZ: $block_size" >>$index_file
        echo -ne "$filter,$nkeys,$nqrys,$keylen,$kdist,$qdist,$minrange,$maxrange,$pqratio,$pnratio,$membudg,$block_size," >>$res_csv
    elif [ $filter = "Oasis" ]; then
        config+="-${membudg}-${block_size}"
        echo -e "BPK: $membudg; BlockSZ: $block_size" >>$index_file
        echo -ne "$filter,$nkeys,$nqrys,$keylen,$kdist,$qdist,$minrange,$maxrange,$pqratio,$pnratio,$membudg,$block_size," >>$res_csv
    fi

    echo -e "\n" >>$index_file
    echo "Running $out_file: " $config

    if [ $filter = "Oasis" ]; then
        echo -e "\tFilter Bits-per-Key:\t$membudg" >>./experiment_result
        echo -e "\tBlock Size:\t$block_size" >>./experiment_result
        echo -e "### END EXPERIMENT DESCRIPTION ###\n\n" >>./experiment_result
        $EXP_BIN "$res_csv" "$filter" "$membudg" "$block_size" >>./experiment_result

    elif [ $filter = "OasisPlus" ]; then
        echo -e "\tFilter Bits-per-Key:\t$membudg" >>./experiment_result
        echo -e "\tElement Per Block:\t$block_size" >>./experiment_result
        echo -e "### END EXPERIMENT DESCRIPTION ###\n\n" >>./experiment_result
        $EXP_BIN "$res_csv" "$filter" "$membudg" "$block_size" "$maxrange" >>./experiment_result
    fi

    echo $out_file " - Complete "
    cp ./experiment_result "$res_folder/$filter/$out_file.txt"
}

# Oasis Experiments
filter="Oasis"
file_cnt=1
mkdir "$res_folder/$filter"
for nkeys in "${nkeys_arr[@]}"; do
    for keylen in "${keylen_arr[@]}"; do
        for nqrys in "${nqrys_arr[@]}"; do
            for i in "${!kdist_arr[@]}"; do
                for j in "${!minrange_arr[@]}"; do
                    for pnratio in "${pnratio_arr[@]}"; do
                        kdist="${kdist_arr[$i]}"
                        qdist="${qdist_arr[$i]}"
                        minrange="${minrange_arr[$j]}"
                        maxrange="${maxrange_arr[$j]}"
                        pqratio="${pqratio_arr[$j]}"

                        for membudg in "${membudg_arr[@]}"; do
                            for block_size in "${block_sizes[@]}"; do
                                printf -v filecnt "%03d" $file_cnt
                                expdir=$(mktemp -d $SCRATCH_DIR/$filter.XXXXXXX)
                                experiment
                                rm -rf $expdir
                                ((file_cnt++))
                            done
                        done
                    done
                done
            done
        done
    done
done


# OasisPlus Experiments
filter="OasisPlus"
file_cnt=1
mkdir "$res_folder/$filter"
for nkeys in "${nkeys_arr[@]}"; do
    for keylen in "${keylen_arr[@]}"; do
        for nqrys in "${nqrys_arr[@]}"; do
            for i in "${!kdist_arr[@]}"; do
                for j in "${!minrange_arr[@]}"; do
                    for pnratio in "${pnratio_arr[@]}"; do
                        kdist="${kdist_arr[$i]}"
                        qdist="${qdist_arr[$i]}"
                        minrange="${minrange_arr[$j]}"
                        maxrange="${maxrange_arr[$j]}"
                        pqratio="${pqratio_arr[$j]}"

                        for membudg in "${membudg_arr[@]}"; do
                            for block_size in "${block_sizes[@]}"; do
                                printf -v filecnt "%03d" $file_cnt
                                expdir=$(mktemp -d $SCRATCH_DIR/$filter.XXXXXXX)
                                experiment
                                rm -rf $expdir
                                ((file_cnt++))
                            done
                        done
                    done
                done
            done
        done
    done
done

rm -rf $SCRATCH_DIR
