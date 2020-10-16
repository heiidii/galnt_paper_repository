import os
import sys
i0=1

def getmappingfromlines(lines,option="quiet"):
    line0=lines[i0]
    headersplit = line0.split()
    mapheader = dict()
    i=0
    for entry in headersplit:
                mapheader[entry] = i
                i+=1
                if option=='verbose':
                  print(i,entry)
    return mapheader

def getfiltereddatafromlines(lines,fields_local=None,dictofarrays=None,nodelete=True):
  mapheader = dict()
  mapheader = getmappingfromlines(lines,"quiet")
  dictofarrays = dict()
  fields_avail = {}
  delete_fields =[]
  max_index = 1
  for field in fields_local:
                dictofarrays[field]=[]
                minfield = fields_local[field][0]
                maxfield = fields_local[field][1]
                if field in mapheader:
                        fields_local[field][2] = mapheader[field]
                        if fields_local[field][2] > max_index:
                                max_index = fields_local[field][2]
                        fields_avail[field] = fields_local[field]
                else:
                        delete_fields.append(field)
  
  dictofarrays["description"]=[]  
  for line in lines[i0+1:]:
                data = line.split()
                if line.find("total_score") != -1: continue
                setTrue = False
                addstring = []
                if len(data) < max_index: continue
                tempofarrays = dict()
                if "description" in mapheader:
                   tempofarrays["description"]=data[mapheader["description"]]
                for field in fields_avail:
                        minfield = fields_avail[field][0]
                        maxfield = fields_avail[field][1]
                        index = fields_avail[field][2]
                        val = float(data[index])
                        if val<maxfield and val>=minfield:
                                setTrue = True
                                addstring.append(str(val))
                                tempofarrays[field] = val
        
                        else:
                                setTrue = False
                                break
                if setTrue:
                        for field in fields_avail:
                                dictofarrays[field].append(tempofarrays[field])
                        if "description" in tempofarrays:
                                dictofarrays["description"].append(tempofarrays["description"])
  return dictofarrays, fields_avail

