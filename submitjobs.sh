year=2017
tag=V9_forAN
#_noOrthogonal

path=/eos/cms/store/group/phys_exotica/monoHiggs/monoHbb/skimmedFiles/merged_bkgrootfiles_splitted/
count=1
python -u SubmitJobs.py processSample.py $year $tag $path $count 
count=2
year=2018
path=/eos/cms/store/group/phys_exotica/monoHiggs/monoHbb/2018_skimmedFiles/merged_bkgrootfiles_splitted/
#path=/eos/cms/store/group/phys_exotica/monoHiggs/monoHbb/2018_skimmedFiles/v12.09_addedPuppiMET_v2_merged/
python -u SubmitJobs.py processSample.py $year $tag $path $count 

