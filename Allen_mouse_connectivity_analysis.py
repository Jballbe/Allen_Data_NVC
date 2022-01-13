from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
mcc = MouseConnectivityCache()

# open up a list of all of the experiments
all_experiments = mcc.get_experiments(dataframe=True) # a dataframe
print(all_experiments)
#take a look at what we know about an experiment with a primary motor injection
print(all_experiments.loc[122642490])
