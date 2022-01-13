from allensdk.core.brain_observatory_cache import BrainObservatoryCache

'''
Demo 
'''

boc = BrainObservatoryCache() #the boc functions : https://alleninstitute.github.io/AllenSDK/allensdk.core.brain_observatory_cache.html

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



