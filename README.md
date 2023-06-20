ValidateXL
==========

![Image of ValidateXL flow chart](./media/validateXL_process_diagram.png)

Validate_xl.py is a module for working with Crosslinking Mass Spectrometry 
data that has been analysed by xQuest v2.1.1. xQuest searches LC-MS/MS data
for lysine specific cross-linked peptides and the cross-linking site. 

ValidateXL assess spectra quality of a spectra-crosslink match to allow 
automated validation of crosslinks that pass a quality control test as 
described in the figure above. To be considered for automatic validation 
spectra must have at least 30% matched ions on both the alpha and beta
peptides in the crosslink. There must be at least one matched ion containing
the crosslinker on both the alpha and beta peptides.

ValidateXL uses the xQuest result file; merged_xquest.xml which is generated
and stored in the xQuest result file after a completed search of LC-MS/MS
data. The fasta file used in the original xQuest search is also necessary.

ValidateXL makes use of click for Command Line Inputs.

Inputs: 
    An input_dir containing the following subdirectories:
    'xq_xml' which contains the merged_xquest.xml
    'fasta' which contains the .fasta file used in the original search
    An output_dir where the resulting CSV files will be generated.

Outputs: 3 CSV files for each result file analysed labelled as follows:
    * identifier_validated_results.csv
    These results have sufficent crosslink peaks and sequence coverage to
    be defined as genuine crosslinks.
    * identifier_rejected_results.csv
    These results have poor annotated sequence coverage and should not be
    regarded as crosslinks.
    * identifier_manual_results.csv
    These results would benefit from manual validation of the spectral 
    quality to confirm the presence of a crosslink.

Validate_xl.py currently works on linux based operating systems

To execute enter the location of the input and output directories into the 
command line. as below:

`python validatexl --input_dir=your_input_dir --output_dir=your_output_dir`


===========================================================================
Author: Juliette M.B James 21/11/2017
This code is offered under a GNU GPLv3 License. For more details on 
licensing terms please see: https://choosealicense.com/licenses/gpl-3.0/
