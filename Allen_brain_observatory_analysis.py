from allensdk.core.brain_observatory_cache import BrainObservatoryCache
import allensdk.brain_observatory.stimulus_info as stim_info

'''
Demo 
'''

boc = BrainObservatoryCache() #the boc functions : https://alleninstitute.github.io/AllenSDK/allensdk.core.brain_observatory_cache.html

#Experiment containers data
targeted_structures = boc.get_all_targeted_structures()
print("all targeted structures: " + str(targeted_structures))
cre_lines=boc.get_all_cre_lines()
print("all cre lines: "+str(cre_lines))
imaging_depths=boc.get_all_imaging_depths()
print("all image depths: "+str(imaging_depths))
reporter_lines=boc.get_all_reporter_lines()
print("all reporter lines: "+str(reporter_lines))
session_types=boc.get_all_session_types()
print("all session types: "+str(session_types))
all_stimuli=boc.get_all_stimuli()
print("all stimuli: "+str(all_stimuli))


# Download experiment containers for Cux2 experiments
cux2_ecs = boc.get_experiment_containers(cre_lines=['Cux2-CreERT2'])
print("Cux2 experiments: %d\n" % len(cux2_ecs))
print("Example experiment container record:")
print(cux2_ecs[0])

##Download Experiments for a Container
cux2_ec_id = cux2_ecs[0]['id'] #the experiment containers of the Cux2 experiments
exps = boc.get_ophys_experiments(experiment_container_ids=[cux2_ec_id]) #store in a electrophysiological data file (nwb file) the data of an ophys experiment (has to mention the ophys_experiment_id
print("Experiments for experiment_container_id %d: %d\n" % (cux2_ec_id, len(exps)))
print(exps)

# pick one of the cux2 experiment containers
cux2_ec_id = cux2_ecs[-1]['id']

# Find the experiment with the static static gratings stimulus
exp = boc.get_ophys_experiments(experiment_container_ids=[cux2_ec_id],
                                stimuli=[stim_info.STATIC_GRATINGS])[0]
print("Experiment with static gratings:")
print(exp)

#data_set = boc.get_ophys_experiment_data(exp['id'])
#got OSError: Unable to open file (truncated file: eof = 10190123, sblock->base_addr = 0, stored_eof = 197846895) : may be because the .h5 is partially downloaded or corrupted

# print out the metadata available in the NWB file
#print(data_set.get_metadata())

##Find Cells of interests ; find experiments for Cells ; Download Experiment Data for a Cell ; Fluorescence Traces ; ROI Masks & analysis ; Neuropil Correction ; ....
#cells = boc.get_cell_specimens() #very long bc there are 63251 cells but works in a jupyter notebook 