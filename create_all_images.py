import requests
import json

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO
import base64

import os

from functools import reduce

server_url = "http://localhost:2000/grassroots/public_backend"
##server_url = "http://localhost:80/grassroots/public_backend"

def lookup_keys(dictionary, keys, default=None):
     return reduce(lambda d, key: d.get(key, default) if isinstance(d, dict) else default, keys.split("."), dictionary)

###################################################################
def searchPhenotypeTrait(listPheno, value):

    name = listPheno[value]['definition']['trait']['so:name']

    return name

###################################################################
def searchPhenotypeUnit(listPheno, value):

    name = listPheno[value]['definition']['unit']['so:name']

    return name

###################################################################
def search_phenotype(list_observations, value):

    found = False
    for i in range(len(list_observations)):

        dic            = list_observations[i]
        phenotype_name = lookup_keys(dic, 'phenotype.variable')
        if  (phenotype_name == value ):
              return True
              break

    return found

############################################################################
def search_phenotype_index(list_observations, value):

    for i in range(len(list_observations)):

        dic            = list_observations[i]
        phenotype_name = lookup_keys(dic, 'phenotype.variable')
        if  (phenotype_name == value ):
              return i

####################################################################
def dict_phenotypes(pheno, plots):

    names = []
    traits = []
    for key in pheno:
        #print("-->", key)
         

        names.append(key)
        traits.append(pheno[key]['definition']['trait']['so:name'])

    phenoDict = dict(zip(names, traits))

    for j in range(len(plots)):
        if ( 'discard' in plots[j]['rows'][0] ):
            pass
        if ('observations' in plots[j]['rows'][0]):
            for k in range(len(plots[j]['rows'][0]['observations'])):
                if ('raw_value' in plots[j]['rows'][0]['observations'][k]):
                    rawValue = plots[j]['rows'][0]['observations'][k]['raw_value']
                if ('corrected_value' in plots[j]['rows'][0]['observations'][k]):
                    rawValue = plots[j]['rows'][0]['observations'][k]['corrected_value']
                if ( type(rawValue) == str):
                    name = plots[j]['rows'][0]['observations'][k]['phenotype']['variable']
                    if( name in phenoDict.keys() ):
                        print("Remove:", phenoDict[name])
                        del phenoDict[name]
    
    return phenoDict    


####################################################################
def get_plot(id):
    plot_request = {
            "services": [{
                "so:name": "Search Field Trials",
                "start_service": True,
                "parameter_set": {
                    "level": "advanced",
                    "parameters": [{
                        "param": "ST Id",
                        "current_value": id
                    }, {
                        "param": "Get all Plots for Study",
                        "current_value": True
                    }, {
                        "param": "ST Search Studies",
                        "current_value": True
                    }]
                }
            }]
        }
    res = requests.post(server_url, data=json.dumps(plot_request))
    return json.dumps(res.json())

####################################################################
def get_plots():
    plot_request = {
            "services": [{
                "so:name": "Search Field Trials",
                "start_service": True,
                "parameter_set": {
                    "level": "advanced",
                    "parameters": [{
                        "param": "ST Id",
                        "current_value": ""
                    }, {
                        "param": "Get all Plots for Study",
                        "current_value": True
                    }, {
                        "param": "ST Search Studies",
                        "current_value": True
                    }]
                }
            }]
        }
    res = requests.post(server_url, data=json.dumps(plot_request))
    return json.dumps(res.json())



####################################################################
def get_all_fieldtrials():
    list_all_ft_request = {
        "services": [
            {
                "so:name": "Search Field Trials",
                "start_service": True,
                "parameter_set": {
                    "level": "simple",
                    "parameters": [
                        {
                            "param": "FT Keyword Search",
                            "current_value": ""
                        },

                        {
                            "param": "FT Study Facet",
                            "current_value": True
                        },
                        {
                            "param": "FT Results Page Number",
                            "current_value": 0
                        },
                        {
                            "param": "FT Results Page Size",
                            "current_value": 500
                        }
                    ]
                }
            }
        ]
    }
    res = requests.post(server_url, data=json.dumps(list_all_ft_request))
    return json.dumps(res.json())


##############################################################################
'''
Create numpy arrays for plotly script. Matrix of raw values and matrix of accession 
'''
def numpy_data(json, pheno, current_name, total_rows, total_columns):

    traitName = searchPhenotypeTrait(pheno, current_name)
    unit      = searchPhenotypeUnit( pheno, current_name)

    dtID= np.dtype(('U', 4))

    row_raw   = np.array([])
    matrix    = np.array([])
    row_acc   = np.array([])
    accession = np.array([])
    plotsIds  = np.array([], dtype=dtID)  #format of strings

    matrices = []

    num_columns = 1
    row    = 1
    column = 1
    #loop throght observations in the same fashion as in old JS code. 
    for j in range(len(json)):
        if ( int( json[j]['row_index'] ) == row ):
            if  (int( json[j]['column_index'] ) == column):
               if column > num_columns:
                   num_columns = column

               if   ( 'discard' in json[j]['rows'][0] ):
                    row_raw  = np.append(row_raw, np.nan )  # use NaN for discarded plots
                    row_acc  = np.append(row_acc, np.nan )  
                    plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
               elif ( 'blank' in json[j]['rows'][0] ):
                    row_raw  = np.append(row_raw, np.nan )  # use NaN for discarded plots
                    row_acc  = np.append(row_acc, np.nan )  
                    plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
      
               elif ( 'observations' in json[j]['rows'][0] ):
                    if( search_phenotype(json[j]['rows'][0]['observations'], current_name) ):
                        indexCurrentPhenotype = search_phenotype_index (json[j]['rows'][0]['observations'], current_name)
                        if ('raw_value' in json[j]['rows'][0]['observations'][indexCurrentPhenotype]):
                            rawValue = json[j]['rows'][0]['observations'][indexCurrentPhenotype]['raw_value']
                        if ('corrected_value' in json[j]['rows'][0]['observations'][indexCurrentPhenotype]):    
                            rawValue = json[j]['rows'][0]['observations'][indexCurrentPhenotype]['corrected_value']
                        row_raw  = np.append(row_raw, rawValue) 
                        row_acc  = np.append(row_acc, json[j]['rows'][0]['material']['accession']) 
                        plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
                    else:
                        row_raw  = np.append(row_raw, np.inf )  # use infinity for N/A data
                        row_acc  = np.append(row_acc, json[j]['rows'][0]['material']['accession'])  
                        plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
               else:
                    if ( 'rows' in json[j] ):
                        row_raw  = np.append(row_raw, np.inf )  # use infinity for N/A data
                        row_acc  = np.append(row_acc, json[j]['rows'][0]['material']['accession'])  
                        plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
         
  
               column+=1
               columns = json[j]['column_index']#

        elif ( int( json[j]['row_index'] ) > row  ):
            if column > num_columns:
                   num_columns = column

            if   ( 'discard' in json[j]['rows'][0] ):
                    row_raw  = np.append(row_raw, np.nan )  
                    row_acc  = np.append(row_acc, np.nan )  
                    plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
            elif   ( 'blank' in json[j]['rows'][0] ):
                    row_raw  = np.append(row_raw, np.nan )  
                    row_acc  = np.append(row_acc, np.nan )  
                    plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
        
            elif ( 'observations' in json[j]['rows'][0] ):
                    if( search_phenotype(json[j]['rows'][0]['observations'], current_name) ):
                        indexCurrentPhenotype = search_phenotype_index (json[j]['rows'][0]['observations'], current_name)
                        if ('raw_value' in json[j]['rows'][0]['observations'][indexCurrentPhenotype]):
                            rawValue = json[j]['rows'][0]['observations'][indexCurrentPhenotype]['raw_value']
                        if ('corrected_value' in json[j]['rows'][0]['observations'][indexCurrentPhenotype]):    
                            rawValue = json[j]['rows'][0]['observations'][indexCurrentPhenotype]['corrected_value']
                        row_raw  = np.append(row_raw, rawValue) 
                        row_acc  = np.append(row_acc, json[j]['rows'][0]['material']['accession']) 
                        plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
                    else:
                        row_raw  = np.append(row_raw, np.inf )
                        row_acc  = np.append(row_acc, json[j]['rows'][0]['material']['accession'])  
                        plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
            else:
                    if ( 'rows' in json[j] ):
                        ##print("rows with no observations------",json[j])
                        row_raw  = np.append(row_raw, np.inf )  # use infinity for N/A data
                        row_acc  = np.append(row_acc, json[j]['rows'][0]['material']['accession'])  
                        plotsIds = np.append(plotsIds, json[j]['rows'][0]['study_index'] )
             

            row+=1
            column=2
            columns = json[j]['column_index']


    #column = columns # use actual number of columns instead of counter
    column = num_columns-1

    if column<columns:
        column=columns
    
    #######print("number of plots and shape check", len(json), row, column, row*(column) )
    if (len(json) != row*column):
        #print("NOT rectangular")
        if(total_columns!=None):
          if(column<total_columns):
             column=total_columns

        # fit odd shape plot into bigger rectangular plot.
        row_raw  = oddShapeValues(   json, row, column, current_name)
        #row_acc  = oddShapeAccession(json, row, column, current_name)
        #plotsIds = oddShapePlotID(   json, row, column, current_name)

    matrices.append(row)
    matrices.append(column)
    matrices.append(row_raw)
    #matrices.append(row_acc)
    matrices.append(traitName)
    matrices.append(unit)
    #matrices.append(plotsIds)
    
    return matrices


#####################################################################################################
def oddShapeValues(arraysJson, rows, columns, phenotype):

    #print("ROWS COLUMNS",rows, columns)    

    matrix = np.zeros((rows,columns))
    matrix[:] = np.nan

    for r in range(len(arraysJson)):
        if  ( 'discard' in arraysJson[r]['rows'][0] ):
            i = int( arraysJson[r]['row_index']    )
            j = int( arraysJson[r]['column_index'] )
            i=i-1
            j=j-1
            matrix[i][j] = np.nan
        elif( 'blank' in arraysJson[r]['rows'][0] ):
            i = int( arraysJson[r]['row_index']    )
            j = int( arraysJson[r]['column_index'] )
            i=i-1
            j=j-1
            #print("row i, column j----------", i, j)
            matrix[i][j] = np.nan
    

        elif ( 'observations' in arraysJson[r]['rows'][0] ):
            i = int( arraysJson[r]['row_index']    )
            j = int( arraysJson[r]['column_index'] )
            i=i-1
            j=j-1
            if( search_phenotype(arraysJson[r]['rows'][0]['observations'], phenotype) ):
                indexCurrentPhenotype = search_phenotype_index (arraysJson[r]['rows'][0]['observations'], phenotype) 
                if ('raw_value' in arraysJson[r]['rows'][0]['observations'][indexCurrentPhenotype]):
                    rawValue = arraysJson[r]['rows'][0]['observations'][indexCurrentPhenotype]['raw_value']
                if ('corrected_value' in arraysJson[r]['rows'][0]['observations'][indexCurrentPhenotype]):
                    rawValue = arraysJson[r]['rows'][0]['observations'][indexCurrentPhenotype]['corrected_value']
                matrix[i][j] = rawValue
            else:
                matrix[i][j] = np.inf
        
        else:
            if('rows' in arraysJson[r]):        #rows field exists but it has no observations!
               i = int( arraysJson[r]['row_index']    )
               j = int( arraysJson[r]['column_index'] )
               i=i-1
               j=j-1
               matrix[i][j] = np.inf  # consider it N/A instead as default discarded (nan)        

    #matrix = np.flipud(matrix)
    #print(matrix)
    matrix  = matrix.flatten()

    return matrix
#####################################################################################################
'''
test rendering seaborn image
'''
def seaborn_plot(numpy_matrix, title, unit, uuid, name):

    sns.set(rc={'figure.figsize':(15.5,6.5)})

    numpy_matrix = np.flipud(numpy_matrix)      # To Match order shown originally in JS code
    notAvailable = np.zeros_like(numpy_matrix)
    discarded    =  np.zeros_like(numpy_matrix)
    indexInf     =  np.where(np.isinf(numpy_matrix))
    indexDiscard =  np.where(np.isnan(numpy_matrix))
    notAvailable[indexInf]   = 1
    discarded[indexDiscard]  = 1
    NA        = np.where(notAvailable < 1, np.nan, notAvailable)
    discarded = np.where(   discarded < 1, np.nan, discarded)
    units = 'Units: '+ unit

    numpy_matrix[indexInf] = np.nan # Replace Inf by NaN

    # Reverse Y ticks and start them from 1
    size  = numpy_matrix.shape
    Y     = size[0]
    Yvals = np.arange(0.5, Y+0.5, 1.0)
    Yaxis = np.arange(1,Y+1)
    Yaxis = np.flip(Yaxis)

    X     = size[1]
    Xvals = np.arange(0, X)
    Xaxis = np.arange(1,X+1)


    maxVal = np.nanmax(numpy_matrix)
    minVal = np.nanmin(numpy_matrix)
    #print(minVal)
    colormap  = sns.light_palette("seagreen", as_cmap=True)
    dark      = sns.dark_palette((260, 75, 60), input="husl")
    sns.heatmap(NA, linewidth=0.5,cmap=dark, cbar=False )

    g = sns.heatmap(numpy_matrix,  vmax=maxVal, vmin=minVal,linewidth=0.5,cmap=colormap, cbar_kws={'label': units}) 
    ##g.set_facecolor('xkcd:black')

    g.patch.set_facecolor('white')
    g.patch.set_edgecolor('black')
    g.patch.set_hatch('xx')


    g.set_xlabel("Columns", fontsize = 14)
    g.set_ylabel("Rows", fontsize = 14)
    g.set_title(title, fontsize = 20)
    g.set_yticks( Yvals )
    g.set_yticklabels(Yaxis, size=10)
    g.tick_params(    axis='y', rotation=0)
    g.set_xticks(Xvals)
    g.set_xticklabels(Xaxis, size=10)
   

    fig = g.get_figure()
    mem = BytesIO()
    fig.savefig(mem, format='png')
    script_dir  = os.path.dirname(__file__)
    results_dir = os.path.join(script_dir, uuid)
    #print ("PATHS", results_dir)

    #fig.savefig('heatnew.png')
    name      = name + '.png'
    imageName = name.replace("/","-")
    ##############print ("file name", imageName)
    fig.savefig(os.path.join (results_dir , imageName))
    fig.clf()
    mem.seek(0)
    image = base64.b64encode(mem.getvalue()) # load the bytes in context as base64
    image = image.decode('utf8')
    mem.close()
     
    return image

################################
################################

id='6082cb1b02700f28764e71d6'
study      = get_plot(id)
study_json = json.loads(study)
plot_data  = study_json['results'][0]['results'][0]['data']['plots']

phenotypes = study_json['results'][0]['results'][0]['data']['phenotypes']

all_studies  = get_all_fieldtrials()
all_studies = json.loads(all_studies)

#print("results", all_studies)
print("TOTAL NUMBER", len(all_studies['results'][0]['results']))
#print("results", all_studies['results'][0]['results'][0]['data'])
#print("---", all_studies['results'][0]['results'][116]['data'])

studies_ids =[]
for i in range(len(all_studies['results'][0]['results'])):
#for i in range(10):
	uuid 	     = all_studies['results'][0]['results'][i]['data']['_id']['$oid']

	#single_study = get_plot(uuid) 
	#single_study = json.loads(single_study)
	
	if 'phenotypes' in all_studies['results'][0]['results'][i]['data']:
		studies_ids.append(uuid)		
	

#studies_ids.remove('6054df3a02700f065c6fcc58') #DFW TKNIL Set 1 JIC-Morley, Harvest 2016. add exception

##studies_ids.remove('619e0cf787a2793484741458') #   DFW Zinc NAM RRes, Harvest 2019
##studies_ids.remove('619e0d970598733f093a0085') #   DFW Zinc NAM Brooms Barn, Harvest 2019 ERROR due rows: [ ]


studies_ids.remove('619e159b87a279348474145b') # DFW Academic Toolkit RRes, Harvest 2021   
studies_ids.remove('6225dfde93b7641e4b5acb85') #  NIAB CSSL AB Glasshouse exp 

##studies_ids.remove('62cc45c34fe23f46bd1aeafa') #  DFW-BTK-H2019-Limagrain
##studies_ids.remove('62cda1e5fc5013010a2bc7ee') #  DFW-BTK-H2019-KWS

# L studies_ids.remove('5dd8009ade68e75a927a8274') #1st vs 3rd wheat take-all resistance trial ERROR rows: [ ]


##studies_ids =[]
##studies_ids.append('603e3e9502700f7faf25dfb4')

#studies_ids.append('62c4380ab9d8cf1e3c12aebe') DFW-BTK-H2019-BASF


#for j in range(5):
for j in range(len(studies_ids)):
    study      = get_plot(studies_ids[j])
    study_json = json.loads(study)

    plot_data  = study_json['results'][0]['results'][0]['data']['plots']
    name       = study_json['results'][0]['results'][0]['data']['so:name']
    print("---------study :", name)
    print("study id :",studies_ids[j])

    phenotypes = study_json['results'][0]['results'][0]['data']['phenotypes']
    dictTraits = dict_phenotypes(phenotypes, plot_data)  #create dictionary of phenotypes.
    default    = list(dictTraits.keys())[0]              #use first one as default. 

    total_rows = study_json['results'][0]['results'][0]['data']['num_rows']
    total_columns = study_json['results'][0]['results'][0]['data']['num_columns']

    phenoList   = list(dictTraits.keys()) 
	
    print(name)
    print("number of observations:", len(phenoList))

    dirName = studies_ids[j]
    if not os.path.exists(dirName):
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ")
    else:    
        print("Directory " , dirName ,  " already exists")

    for k in range(len(phenoList)):
        selected_phenotype = phenoList[k]
        ###print(selected_phenotype)
        matrices  = numpy_data(plot_data, phenotypes, selected_phenotype, total_rows, total_columns)

        row     = matrices[0]
        column  = matrices[1]
        row_raw = matrices[2]
        #row_acc = matrices[3]
        traitName = matrices[3]
        units     = matrices[4]
        
        matrix   = row_raw.reshape(row,column)
        image = seaborn_plot(matrix, traitName, units, dirName, selected_phenotype)


