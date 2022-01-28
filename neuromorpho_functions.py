import neuromorpholib
from neuromorpholib import neuromorpho
from optparse import OptionParser
import inspect
nmo = neuromorpho.NeuroMorpho()
print(inspect.getmembers(nmo,predicate=inspect.ismethod))
#mouse_neurons_visp= nmo.search({"species":"mouse","brain_region":"primary visual"}) #34 pages for 1714 cells in visp ; returns a list of dico
#nmo.seach keys : ['neuron_id', 'neuron_name', 'archive', 'age_scale', 'gender', 'reference_pmid', 'reference_doi', 'note', 'age_classification', 'brain_region', 'cell_type', 'species', 'strain', 'scientific_name', 'stain', 'experiment_condition', 'protocol', 'slicing_direction', 'reconstruction_software', 'objective_type', 'original_format', 'domain', 'attributes', 'magnification', 'upload_date', 'deposition_date', 'shrinkage_reported', 'shrinkage_corrected', 'reported_value', 'reported_xy', 'reported_z', 'corrected_value', 'corrected_xy', 'corrected_z', 'soma_surface', 'surface', 'volume', 'slicing_thickness', 'min_age', 'max_age', 'min_weight', 'max_weight', 'png_url', 'Physical_Integrity']
#print(mouse_neurons_visp[0])



#key of nmo list : 'neuron_id , 'neuron_name' of the file, 'archive' gives the author (Allen Cell Types, etc.), 'note' , ' age_scale' , 'gender', 'age_classification' adult, 'brain_region'
#'brain_region' : ['neocortex', 'occipital', 'primary visual', 'layer 4'], 'cell_type': ['interneuron', 'Aspiny'], 'species': 'mouse', 'strain': 'Chat-IRES-Cre-neo'
#scientific_name' , 'stain' , 'experiment_condition' , 'protocol' , 'slicing_direction' , 'reconstruction_software' , ' objective_type' , 'original_format', 'domain': 'Dendrites, Soma, Axon',
# 'attributes', 'magnification', 'upload_date', 'deposition_date', 'shrinkage_reported', 'shrinkage_corrected','reported_value', 'reported_xy', 'reported_z', 'corrected_value', 'corrected_xy',
# 'corrected_z', 'soma_surface', 'surface', 'volume', 'slicing_thickness', 'min_age', 'max_age', 'min_weight', 'max_weight', 'png_url' , 'reference_pmid', 'reference_doi', '
# physical_Integrity', '_links', 'measurements', 'persistence_vector' .


def neuromorpho_layer (structure_name) : #primary visual
    '''
    For a given structure it returns the index of all its neurons according to the different layers

    Parameters
    ----------
    structure_name : str
        name of the structure : 'primary visual'

    Returns
    -------
    l1, l2_3, l3_4, l4, l5, l6, l6a, l6b, no_layer : list
        List of index of the neurons associated to the different layers

    '''

    visp = nmo.search({"species": "mouse", "brain_region": structure_name})
    l1 = []
    l2_3 = []
    l4 = []
    l5 = []
    l6 = []
    l6a = []
    l6b = []
    no_layer = []
    l3_4 = []
    other = []
    region_other = []

    index = 0
    for i in visp:

        if i["brain_region"] == ['neocortex', 'occipital', 'primary visual', 'layer 1']:
            l1.append(index)
        else:
            if i["brain_region"] == ['neocortex', 'occipital', 'primary visual', 'layer 2-3'] or i["brain_region"] == [
                'neocortex', 'occipital', 'primary visual', 'layer 2-3', 'right'] or i["brain_region"] == ['neocortex',
                                                                                                           'occipital',
                                                                                                           'primary visual',
                                                                                                           'layer 2'] or \
                    i["brain_region"] == ['neocortex', 'occipital', 'primary visual', 'layer 3'] or i[
                'brain_region'] == ['neocortex', 'occipital', 'primary visual', 'layer 2-3 Binocular region']:
                l2_3.append(index)
            else:
                if i["brain_region"] == ['neocortex', 'occipital', 'primary visual', 'layer 4'] or i[
                    "brain_region"] == ['neocortex', 'occipital', 'primary visual', 'layer 4', 'left'] or i[
                    "brain_region"] == ['neocortex', 'occipital', 'primary visual', 'layer 4', 'right']:
                    l4.append(index)
                else:
                    if i["brain_region"] == ['neocortex', 'occipital', 'primary visual', 'layer 5'] or i[
                        "brain_region"] == ['neocortex', 'occipital', 'primary visual', 'layer 5', 'left'] or i[
                        "brain_region"] == ['neocortex', 'occipital', 'primary visual', 'layer 5', 'right'] or i[
                        "brain_region"] == ['neocortex', 'occipital', 'primary visual', 'layer 5, left'] or i[
                        'brain_region'] == ['neocortex', 'occipital', 'primary visual', 'layer 5 Binocular region'] or \
                            i['brain_region'] == ['neocortex', 'occipital', 'primary visual', 'layer 5 left']:
                        l5.append(index)
                    else:
                        if i["brain_region"] == ['neocortex', 'occipital', 'primary visual', 'layer 6',
                                                 'superficial'] or i["brain_region"] == ['neocortex', 'occipital',
                                                                                         'primary visual',
                                                                                         'layer 6a'] or i[
                            "brain_region"] == ['neocortex', 'occipital', 'primary visual', 'layer 6a', 'left']:
                            l6a.append(index)
                        else:
                            if i["brain_region"] == ['neocortex', 'occipital', 'primary visual', 'layer 6', 'deep'] or \
                                    i["brain_region"] == ['neocortex', 'occipital', 'primary visual', 'layer 6b']:
                                l6b.append(index)
                            else:
                                if i["brain_region"] == ['neocortex', 'occipital', 'primary visual']:
                                    no_layer.append(index)
                                else:
                                    if i["brain_region"] == ['neocortex', 'occipital', 'primary visual', 'layer 3-4']:
                                        l3_4.append(index)
                                    else:
                                        if i['brain_region'] == ['neocortex', 'occipital', 'primary visual', 'layer 6']:
                                            l6.append(index)
                                        else:
                                            other.append(index)
                                            region_other.append(i["brain_region"])
        index = index + 1
    return(l1,l2_3,l3_4,l4,l5,l6,l6a,l6b,no_layer)

a,b,c,d,e,f,g,h,j =neuromorpho_layer('primary visual')
print(len(a),len(b),len(c),len(d),len(e),len(f),len(g),len(h),len(j))


#écrire fonction qui retourne la liste d'id à télécharger pour chaque layer pour une brain region mentionnée (primary visual), si exci/inhi mentionné : faire un dico où pour chaque layer on a
#la liste d'id pour exci/inhi(indéterminé

#exci = ['principal cell', 'Spiny'] , ['principal cell', 'pyramidal', 'projection'] , ['principal cell', 'pyramidal'] , ['principal cell', 'pyramidal', 'Excitatory']
# ['principal cell', 'Spiny', 'Sodium Channel', 'Non Voltage Gated 1 Alpha Subunit-positive'] , ['principal cell', 'Spiny', 'Plasma Retinol-Binding Protein-Positive'],
# ['principal cell', 'Spiny', 'Steroid Hormone Receptor Ad4BP-positive'] , ['principal cell', 'Spiny', 'Connective Tissue Growth Factor positive'] ,
# ['principal cell', 'Spiny', 'Homeobox Protein Cux-2-positive'] , ['principal cell', 'pyramidal', 'glutamatergic'] , ['pyramidal', 'principal cell'] , ['Spiny', 'principal cell']
#['Non Voltage Gated 1 Alpha Subunit-positive', 'Sodium Channel', 'Spiny', 'principal cell'] , 'Retinoid-Related Orphan Receptor-Beta-positive', 'Spiny', 'principal cell'] ,
#['Plasma Retinol-Binding Protein-Positive', 'Spiny', 'principal cell'] , ['Neurotensin Receptor 1-positive', 'Spiny', 'principal cell'] ,
# ['Steroid Hormone Receptor Ad4BP-positive', 'Spiny', 'principal cell'] , ['principal cell', 'pyramidal', 'Calcium/calmodulin-dependent protein kinase II (CaMKII)-positive'] ,
#['principal cell', 'pyramidal', 'projection', 'Intratelencephalic (IT)'] , ['principal cell', 'pyramidal', 'Thick-tufted', 'lateral posterior thalamic nucleus-projecting'] ,

#inhi = ['interneuron', 'Aspiny'] or ['interneuron', 'Parvalbumin (PV)-positive'] or ['interneuron', 'Aspiny', 'Somatostatin (SOM)-positive'] ,
# ['interneuron', 'Aspiny', 'Plasma Retinol-Binding Protein-Positive'] or ['interneuron', 'Aspiny', 'Parvalbumin (PV)-positive'] ,['interneuron', 'Aspiny', 'Steroid Hormone Receptor Ad4BP-positive'],
#['interneuron', 'Aspiny', 'serotonergic'] or ['interneuron', 'Aspiny', 'Serotonin receptor type 3A(5HT3)-positive'] or ['interneuron', 'deep projecting cell'] ,
# ['interneuron', 'horizontally elongated'] or ['interneuron', 'shrub cell'] or ['interneuron', 'basket'] or ['interneuron', 'neurogliaform'] or ['interneuron', 'Martinotti'] ,
# ['interneuron', 'Chandelier'] or ['interneuron', 'double bouquet'] or ['interneuron', 'bipolar'] or ['interneuron', 'bitufted'] or ['interneuron', 'neurogliaform', 'elongated'] ,
# ['interneuron', 'single bouquet'] or ['interneuron', 'basket', 'Type 1', 'Parvalbumin (PV)-positive'] or ['interneuron', 'basket', 'Type 2', 'Parvalbumin (PV)-positive'] ,
# ['interneuron', 'Martinotti', 'Somatostatin (SOM)-positive'] or ['Local projecting', 'Fast-spiking', 'interneuron'] or ['Translaminar', 'Fast-spiking', 'interneuron'] ,
# ['Somatostatin (SOM)-positive', 'Aspiny', 'interneuron'] or ['Non Voltage Gated 1 Alpha Subunit-positive', 'Sodium Channel', 'Spiny', 'principal cell'] ,
#['Non Voltage Gated 1 Alpha Subunit-positive', 'Sodium Channel', 'Aspiny', 'interneuron'] or ['Retinoid-Related Orphan Receptor-Beta-positive', 'Aspiny', 'interneuron'] ,
#['Parvalbumin (PV)-positive', 'Aspiny', 'interneuron'] or ['Aspiny', 'interneuron'] or ['Serotonin receptor type 3A(5HT3)-positive', 'Aspiny', 'interneuron'] or ['GABAergic', 'Aspiny', 'interneuron'] ,
#['Fast-spiking', 'basket', 'interneuron'] or ['interneuron', 'bipolar', 'Vasoactive Intestinal Peptide (VIP)-positive'] ,
# ['interneuron', 'basket', 'Fast-spiking', 'Parvalbumin (PV)-positive', 'double-bouquet'] or ['interneuron', 'basket', 'Horizontal', 'Fast-spiking', 'Parvalbumin (PV)-positive'] ,
#['interneuron', 'basket', 'Large', 'Fast-spiking', 'Parvalbumin (PV)-positive'] or ['interneuron', 'basket', 'Small', 'Fast-spiking', 'Parvalbumin (PV)-positive'] ,

#tumoral : ['principal cell', 'pyramidal', 'peritumoral'] ,

#sparsely spiny : ['principal cell', 'sparsely spiny'] ,

#indéterminé : ['principal cell', 'projection']

#in L5 exci with ['Somatostatin (SOM)-positive', 'Spiny', 'principal cell']

#exc

def exci_inhi (structure_name, index_list) :
    visp = nmo.search({"species": "mouse", "brain_region": structure_name})
    exci=[]
    inhi=[]
    tumoral =[]
    sparsely=[]
    indeterminate=[]
    others=[]
    index=0
    for i in index_list :

        if visp[i]['cell_type'] == ['principal cell', 'Spiny'] or visp[i]['cell_type']==['principal cell', 'pyramidal', 'projection'] or visp[i]['cell_type']==['principal cell', 'pyramidal'] or visp[i]['cell_type']==['principal cell', 'pyramidal', 'Excitatory'] or visp[i]['cell_type']==['principal cell', 'Spiny', 'Sodium Channel', 'Non Voltage Gated 1 Alpha Subunit-positive'] or visp[i]['cell_type']==['principal cell', 'Spiny', 'Plasma Retinol-Binding Protein-Positive'] or visp[i]['cell_type']==['principal cell', 'Spiny', 'Steroid Hormone Receptor Ad4BP-positive'] or visp[i]['cell_type']==['principal cell', 'Spiny', 'Connective Tissue Growth Factor positive'] or visp[i]['cell_type']==['principal cell', 'Spiny', 'Homeobox Protein Cux-2-positive'] or visp[i]['cell_type']==['principal cell', 'pyramidal', 'glutamatergic'] or visp[i]['cell_type']==['pyramidal', 'principal cell'] or visp[i]['cell_type']==['Spiny', 'principal cell'] or visp[i]['cell_type']==['Non Voltage Gated 1 Alpha Subunit-positive', 'Sodium Channel', 'Spiny', 'principal cell'] or visp[i]['cell_type']==['Retinoid-Related Orphan Receptor-Beta-positive', 'Spiny', 'principal cell'] or visp[i]['cell_type']==['Plasma Retinol-Binding Protein-Positive', 'Spiny', 'principal cell'] or visp[i]['cell_type']==['Neurotensin Receptor 1-positive', 'Spiny', 'principal cell'] or visp[i]['cell_type']==['Steroid Hormone Receptor Ad4BP-positive', 'Spiny', 'principal cell'] or visp[i]['cell_type']==['principal cell', 'pyramidal', 'Calcium/calmodulin-dependent protein kinase II (CaMKII)-positive'] or visp[i]['cell_type']==['principal cell', 'pyramidal', 'projection', 'Intratelencephalic (IT)'] or visp[i]['cell_type']==['principal cell', 'pyramidal', 'Thick-tufted', 'lateral posterior thalamic nucleus-projecting'] or visp[i]['cell_type']==['principal cell', 'Excitatory'] or visp[i]['cell_type']==['Somatostatin (SOM)-positive', 'Spiny', 'principal cell'] or visp[i]['cell_type']==['principal cell', 'pyramidal', 'projection', 'Extratelencephalic (ET)']:
            exci.append(index)
        else :
            if visp[i]['cell_type'] == ['interneuron', 'Aspiny'] or visp[i]['cell_type']==['interneuron', 'Parvalbumin (PV)-positive'] or visp[i]['cell_type']==['interneuron', 'Aspiny', 'Somatostatin (SOM)-positive'] or visp[i]['cell_type']==['interneuron', 'Aspiny', 'Plasma Retinol-Binding Protein-Positive'] or visp[i]['cell_type']==['interneuron', 'Aspiny', 'Parvalbumin (PV)-positive'] or visp[i]['cell_type']==['interneuron', 'Aspiny', 'Steroid Hormone Receptor Ad4BP-positive']or visp[i]['cell_type']==['interneuron', 'Aspiny', 'serotonergic'] or visp[i]['cell_type']==['interneuron', 'Aspiny', 'Serotonin receptor type 3A(5HT3)-positive'] or visp[i]['cell_type']==['interneuron', 'deep projecting cell'] or visp[i]['cell_type']==['interneuron', 'horizontally elongated'] or visp[i]['cell_type']==['interneuron', 'shrub cell'] or visp[i]['cell_type']==['interneuron', 'basket'] or visp[i]['cell_type']==['interneuron', 'neurogliaform'] or visp[i]['cell_type']==['interneuron', 'Martinotti'] or visp[i]['cell_type']==['interneuron', 'Chandelier'] or visp[i]['cell_type']==['interneuron', 'double bouquet'] or visp[i]['cell_type']==['interneuron', 'bipolar'] or visp[i]['cell_type']==['interneuron', 'bitufted'] or visp[i]['cell_type']==['interneuron', 'neurogliaform', 'elongated'] or visp[i]['cell_type']==['interneuron', 'single bouquet'] or visp[i]['cell_type']==['interneuron', 'basket', 'Type 1', 'Parvalbumin (PV)-positive'] or visp[i]['cell_type']==['interneuron', 'basket', 'Type 2', 'Parvalbumin (PV)-positive'] or visp[i]['cell_type']==['interneuron', 'Martinotti', 'Somatostatin (SOM)-positive'] or visp[i]['cell_type']==['Local projecting', 'Fast-spiking', 'interneuron'] or visp[i]['cell_type']==['Translaminar', 'Fast-spiking', 'interneuron'] or visp[i]['cell_type']==['Somatostatin (SOM)-positive', 'Aspiny', 'interneuron'] or visp[i]['cell_type']==['Non Voltage Gated 1 Alpha Subunit-positive', 'Sodium Channel', 'Spiny', 'principal cell'] or visp[i]['cell_type']==['Non Voltage Gated 1 Alpha Subunit-positive', 'Sodium Channel', 'Aspiny', 'interneuron'] or visp[i]['cell_type']==['Retinoid-Related Orphan Receptor-Beta-positive', 'Aspiny', 'interneuron'] or visp[i]['cell_type']==['Parvalbumin (PV)-positive', 'Aspiny', 'interneuron'] or visp[i]['cell_type']==['Aspiny', 'interneuron'] or visp[i]['cell_type']==['Serotonin receptor type 3A(5HT3)-positive', 'Aspiny', 'interneuron'] or visp[i]['cell_type']==['GABAergic', 'Aspiny', 'interneuron'] or visp[i]['cell_type']==['Fast-spiking', 'basket', 'interneuron'] or visp[i]['cell_type']==['interneuron', 'bipolar', 'Vasoactive Intestinal Peptide (VIP)-positive'] or visp[i]['cell_type']==['interneuron', 'basket', 'Fast-spiking', 'Parvalbumin (PV)-positive', 'double-bouquet'] or visp[i]['cell_type']==['interneuron', 'basket', 'Horizontal', 'Fast-spiking', 'Parvalbumin (PV)-positive'] or visp[i]['cell_type']==['interneuron', 'basket', 'Large', 'Fast-spiking', 'Parvalbumin (PV)-positive'] or visp[i]['cell_type']==['interneuron', 'basket', 'Small', 'Fast-spiking', 'Parvalbumin (PV)-positive'] or visp[i]['cell_type']==['principal cell', 'Aspiny'] or visp[i]['cell_type']==['interneuron', 'basket', 'Parvalbumin (PV)-positive']:
                inhi.append(index)
            else :
                if visp[i]['cell_type'] == ['principal cell', 'pyramidal', 'peritumoral'] :
                    tumoral.append(index)
                else :
                    if visp[i]['cell_type'] == ['principal cell', 'sparsely spiny'] :
                        sparsely.append(index)
                    else :
                        if visp[i]['cell_type'] == ['principal cell', 'projection'] :
                            indeterminate.append(index)
                        else :
                            others.append(index)
                            print(visp[i]['cell_type'])
        index = index + 1
    return(exci,inhi,tumoral,sparsely,indeterminate,others)

ex,inh,tum,spars,indeter,others=exci_inhi('primary visual',d)
print(len(ex),len(inh),len(tum),len(spars),len(indeter),len(others))
