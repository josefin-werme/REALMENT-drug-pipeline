- remove / replace readme
- remove the check_imputed / check_pred dirs ? (or keep for completion)
- where is predict_drugs.check_completed.R ?? 
- check the "collated="$grex/all_collated" # all_ind collated into single df (and class N subset)"
--- is this really everyone? if so is it per gene?
- where to put utils scripts? 
- function rerun() in impute_grex.job

# signatures
- add data/signatures/filter_lincs_exemplar_mapped.R
- add exemplar_signatures_mappedIDs.txt ? 







# TODO: 
- check scripts what files are taken from data/signatures
- get LINCS download/process code from local file ~/DrugSignatures/signatureSearchProcessData.R)
  
# (from signatures readme)
-----
# OLD FILES
### lincs2.h5 (OLD)
This contains all ~135k exemplar signatures (see comments in local file ~/DrugSignatures/signatureSearchProcessData.R)

### compoundIDs_filtered.txt (OLD)
First filtering done to remove presumed experimental compounds that only had a 'code' and did not have a different cmap or alias
Was used in drug prediction script to filter the lincs data (resulting in ~60k signatures)

# NEW FILES
### exemplar_signatures_mappedIDs.txt (CURRENT)
Examplar compounds where the IDs have been mapped to drug_names (see ~/DrugSignatures/filter_signatures.R), and all compounds that weren't mappable have been excluded

### filter_lincs_exemplar_mapped.R
Filters lincs2.h5 to signatures in exemplar_signatures_mappedIDs.txt

### lincs2_mapped.dat
Filtered to signatures in exemplar_signatures_mappedIDs.txt

### atc_processed.txt
Formatted atc file in drug environment (drug$atc) produced by read_atc() function in local script "~/DrugSignatures/* validate.R"
-----
