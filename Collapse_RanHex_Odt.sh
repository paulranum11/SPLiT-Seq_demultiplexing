#! /bin/bash

#Example Results Directory
#cp /Users/ranump/Desktop/Charlie_test/SPLiT-Seq_demultiplexing/results-UMI


RanHex_barcodes="./RanHex.txt"
Odt_barcodes="./OligoDT.txt"
OUTPUT_DIR="results"

# Read in the RandomHexamer and OligoDT RT primer barcode sequences
declare -a RanHex_BARCODES=( $(cut -b 1- $RanHex_barcodes) )
declare -a Odt_BARCODES=( $(cut -b 1- $Odt_barcodes) )
declare -a Cells_List=( $(ls $OUTPUT_DIR/) )

ls $OUTPUT_DIR > Cells_List.txt

# Walk through all RanHex and Odt barcodes in pairs
# These pairs reflect the barcode pairs added to the first 48 wells of the ROUND1 plate
counter=0
for barcode1 in "${RanHex_BARCODES[@]}";
    do
    RhxBarcode="$barcode1"
    #echo "RhxBarcode is $RhxBarcode"
    OdtBarcode="${Odt_BARCODES[$counter]}"
    
    counter=$((counter+1))

    #Begin finding and collapsing cells 
    #Find all cells with the OligoDT barcode.
    declare -a OdtBarcode_hits=( $(grep "$OdtBarcode" Cells_List.txt) )
    #Find all cells with the RanHex barcode.
    declare -a RhxBarcode_hits=( $(grep "$RhxBarcode" Cells_List.txt) )
    for Odt_hit in "${OdtBarcode_hits[@]}";
        do
        #echo "OligoDT_hit is $Odt_hit"
        #If the OligoDT barcode is in the 1st position return true
        if [ $(echo $Odt_hit | awk 'BEGIN { FS = "-" } ; { print $1 }') = $OdtBarcode ]
        then
            #echo "Odt hits were found in the list"
            for Rhx_hit in "${RhxBarcode_hits[@]}";
                do
                #echo "Random Hexamer hit is $Rhx_hit"
                if [ $(echo $Rhx_hit | awk 'BEGIN { FS = "-" } ; { print $1 }') = $RhxBarcode ]
                then
                	#At this point we have found cells primed with the pair of Odt and Rhx primers from a single well
                	#Next we need to match these reads on the rest of their barcodes.
                	Rhx_bc23=( $(echo "$Rhx_hit" | awk 'BEGIN { FS = "-" } ; { print $2 "-" $3 }') )
                	Odt_bc23=( $(echo "$Odt_hit" | awk 'BEGIN { FS = "-" } ; { print $2 "-" $3 }') )
                	if [ $Rhx_bc23 = $Odt_bc23 ]
                	then
                        echo "collapsing"
                        echo "$Odt_hit and"
                        echo "$Rhx_hit"
                		cat $OUTPUT_DIR/$Rhx_hit >> $OUTPUT_DIR/$Odt_hit
                	    rm $OUTPUT_DIR/$Rhx_hit
                    fi
                fi
            	done
        fi
		done
   		
	done
echo "Cleaning up"
rm Cells_List.txt
