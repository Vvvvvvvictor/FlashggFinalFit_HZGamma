# Hold defs of writing functions for datacard
import os, sys, re
from commonTools import *
from commonObjects import *

def writePreamble(f,options):
  f.write("CMS HGG Datacard - %s - %s\n"%(options.years,sqrts__)) 
  f.write("Auto-generated by flashggFinalFits/Datacard/makeDatacard.py\n")
  f.write("Run with: combine\n")
  f.write("---------------------------------------------\n")
  f.write("imax *\n")
  f.write("jmax *\n")
  f.write("kmax *\n")
  f.write("---------------------------------------------\n")
  return True

def writeProcesses(f,d,options):
  f.write("\n")
  # If opt.prune then remove all rows from dataFrame with prune=1
  if options.prune: d = d[d['prune']==0]
  # d = Pandas DataFrame
  # Shapes
  # Loop over categories in dataframe
  years = list(set([re.search(r'\d{4}', item).group() for item in options.years.split(",") if re.search(r'\d{4}', item)]))
  for cat in d.cat.unique():
    # Loop over rows for respective category
    for year in years:
      for ir,r in d[d['cat']==cat].iterrows():
        if (not (str(year) in r["year"])) and (not (r["year"] == "merged")): continue
        # Write to datacard
        if r['proc'] == "bkg_mass":
          if options.mergeYears:
            if year == years[-1]:
              index = r['model'].rfind("_13TeV_bkgshape")
              if index != -1:
                model = r['model'][:index] + "_13TeV_bkgshape"
                f.write("shapes      %-55s %-40s %s %s\n"%(r['proc'],r['cat'],r['modelWSFile'],model))
          else:
            index = r['model'].rfind("_13TeV_bkgshape")
            if index != -1:
              model = r['model'][:index] + "_" + year + "_13TeV_bkgshape"
              f.write("shapes      %-55s %-40s %s %s\n"%(r['proc'],r['cat'],r['modelWSFile'],model))
        elif (r['proc'] == "data_obs") and options.mergeYears:
          if year == years[-1]:
            f.write("shapes      %-55s %-40s %s %s\n"%(r['proc'],r['cat'],r['modelWSFile'],r['model']))
        else:
          f.write("shapes      %-55s %-40s %s %s\n"%(r['proc'],r['cat'],r['modelWSFile'],r['model']))

  # Bin, observation and rate lines
  lbreak = '----------------------------------------------------------------------------------------------------------------------------------'
  lbin_cat = '%-30s'%"bin"
  lobs_cat = '%-30s'%"observation"
  lbin_procXcat = '%-30s'%"bin"
  lproc = '%-30s'%"process"
  lprocid = '%-30s'%"process"
  lrate = '%-30s'%"rate"        
  # Loop over categories
  for cat in d.cat.unique():
    lbin_cat += "%-55s "%cat
    lobs_cat += "%-55s "%"-1"
    sigID = 0
    # Loop over rows for respective category
    for ir,r in d[d['cat']==cat].iterrows():
      if r['proc'] == "data_obs": continue
      lbin_procXcat += "%-55s "%cat
      lproc += "%-55s "%r['proc']
      if r['proc'] == "bkg_mass": lprocid += "%-55s "%"1"
      else:
        lprocid += "%-55s "%sigID
        sigID -= 1
      if r['rate'] == 1.0: lrate += "%-55.1f "%r['rate']
      else: lrate += "%-55.7f "%r['rate']
  #Remove final space from lines and add to file
  f.write("\n")
  for l in [lbreak,lbin_cat,lobs_cat,lbreak,lbin_procXcat,lproc,lprocid,lrate,lbreak]: 
    l = l[:-1]
    f.write("%s\n"%l)
    
  f.write("\n")
  return True


def writeSystematic(f,d,s,options,stxsMergeScheme=None,scaleCorrScheme=None):

  # For signal shape systematics add simple line
  if s['name'] == 'energyErrShift':
    print("yes")
  if s['type'] == 'signal_shape':
    stitle = "%s_%s"%(outputWSNuisanceTitle__,s['title'])
    if s['mode'] != 'other': stitle += "_%s"%outputNuisanceExtMap[s['mode']]
    # If not correlated: separate nuisance per year
    if s['mode'] in ['scales','smears']:
      for year in options.years.split(","):
        stitle_y = "%s_%s"%(stitle,year) 
        lsyst = "%-70s  param    %-6s %-6s"%(stitle_y,s['mean'],s['sigma'])
        f.write("%s\n"%lsyst)
    else:
      lsyst = "%-70s  param    %-6s %-6s"%(stitle,s['mean'],s['sigma'])
      f.write("%s\n"%lsyst)
    return True
 
  # Else: for yield variation uncertainties...
  # Remove all rows from dataFrame with prune=1 (includes NoTag)
  if options.prune:
    mask = (d['prune']==0)
    d = d[mask]

  # If theory: loop over tiers else run over once
  tiers = []
  if 'tiers' in s: tiers = s['tiers']
  if(not options.doSTXSMerging)&('mnorm' in tiers): tiers.remove("mnorm")
  if len(tiers)==0: tiers = ['']
  for tier in tiers:
    if tier != '': tierStr = "_%s"%tier
    else: tierStr = ''
    
    # If calculating merged bin: loop over mergings else run over once
    mns = []
    if tier == 'mnorm':
      if options.doSTXSMerging:
        for mergeName in stxsMergeScheme: mns.append(mergeName)
    if len(mns) == 0: mns.append('')
    for mn in mns:
      if mn != '': mergeStr = "_%s"%mn
      else: mergeStr = ''
    
      # Construct syst line/lines if separate by year
      if(s['correlateAcrossYears'] == 1)|(s['correlateAcrossYears'] == -1):
	stitle = "%s%s%s"%(s['title'],mergeStr,tierStr)
	# print("stitle", stitle)
	lsyst = '%-50s  %-10s    '%(stitle,s['prior'])
	# Loop over categories and then iterate over rows in category
	for cat in d.cat.unique():
	  for ir,r in d[d['cat']==cat].iterrows():
	    if r['proc'] == "data_obs": continue
	    # Extract value and add to line (with checks)
	    sval = r["%s%s%s"%(s['name'],mergeStr,tierStr)]
	    # if stitle == "CMS_hgg_SigmaEOverEShift": print("%s%s%s"%(s['name'],mergeStr,tierStr))
	    # if stitle == "CMS_hgg_SigmaEOverEShift": print(sval)
	    lsyst = addSyst(lsyst,sval,stitle,r['proc'],cat,r['numEvents'])
	# Remove final space from line and add to file
	f.write("%s\n"%lsyst[:-1])
        # For uncorrelated scale weights: not for merged bins
        if options.doSTXSScaleCorrelationScheme:
          if(tier!='mnorm')&("scaleWeight" in s['name']):
	    for ps,psProcs in scaleCorrScheme.iteritems():
	      psStr = "_%s"%ps
	      stitle = "%s%s%s"%(s['title'],psStr,tierStr)
	      lsyst = '%-50s  %-10s    '%(stitle,s['prior'])
	      # Loop over categories and then iterate over rows in category
	      for cat in d.cat.unique():
		for ir,r in d[d['cat']==cat].iterrows():
		  if r['proc'] == "data_obs": continue
		  # Remove year+hgg tags from proc
		  p = re.sub("_2016_hgg","",r['proc'])
		  p = re.sub("_2017_hgg","",p)
		  p = re.sub("_2018_hgg","",p)
		  # Add value if in proc in phase space else -
		  print("sname", s['name'], "here", p, psProcs) 
		  if p in psProcs: sval = r["%s%s"%(s['name'],tierStr)]
		  else: sval = '-'
		  lsyst = addSyst(lsyst,sval,stitle,r['proc'],cat,r['numEvents'])
              # Remove final space from line and add to file
              f.write("%s\n"%lsyst[:-1])
      else:
	for year in options.years.split(","):
	  stitle = "%s%s%s_%s"%(s['title'],mergeStr,tierStr,year)
	  sname = "%s%s%s_%s"%(s['name'],mergeStr,tierStr,year)
	  lsyst = '%-50s  %-10s    '%(stitle,s['prior'])
	  # Loop over categories and then iterate over rows in category
	  for cat in d.cat.unique():
	    for ir,r in d[d['cat']==cat].iterrows():
	      if r['proc'] == "data_obs": continue
	      # Extract value and add to line (with checks)
	      sval = r[sname]
	      # print("here", sval, sname)
	      lsyst = addSyst(lsyst,sval,stitle,r['proc'],cat,r['numEvents'])
	  # Remove final space from line and add to file
	  f.write("%s\n"%lsyst[:-1])
  return True
          

def addSyst(l,v,s,p,c,n):
  #l-systematic line, v-value, s-systematic title, p-proc, c-cat, n-numEvents
  minEventNumber = 100
  if type(v) is str: 
    l += "%-15s "%v
    return l
  elif type(v) is list: 
    # Symmetric:
    if len(v) == 1: 
      # Check 1: variation is non-negligible. If not then skip
      if abs(v[0]-1)<0.0005: l += "%-15s "%"-"
      # If number of events smaller than minEventNumber add - to datacard
      elif (n < minEventNumber) and (not ('lumi' in s)) and (not ('scale' in s)): l += "%-15s "%"-"
      # Check 2: variation is not negative. Print message and add - to datacard (cleaned later)
      elif v[0] < 0.: 
        print " --> [WARNING] systematic %s: negative variation for (%s,%s)"%(s,p,c)
        #vstr = "%s"%v[0]
        vstr = "-"
        # l += "%-15s "%v[0] # After discussion on the 10.04.24
        l += "%-15s "%vstr
      else:
        l += "%-15.3f "%v[0]
    # Anti-symmetirc
    if len(v) == 2:
      # Check 1: variation is non-negligible. If not then skip
      if(abs(v[0]-1)<0.0005)&(abs(v[1]-1)<0.0005): l += "%-15s "%"-"
      # If number of events smaller than minEventNumber add - to datacard
      elif (n < minEventNumber) and (not ('lumi' in s)) and (not ('scale' in s)): l += "%-15s "%"-"
      # Check 2: neither variation is negative. Print message and add - to datacard (cleaned later)
      elif(v[0]<0.)|(v[1]<0.):
        print " --> [WARNING] systematic %s: negative variation for (%s,%s)"%(s,p,c)
        #vstr = "%.3f/%.3f"%(v[0],v[1])
        vstr = "-"
        l += "%-15s "%vstr
      # Check 3: effect is approximately symmetric: then just add single up variation
      elif( abs((v[0]*v[1])-1)<0.0005 ): l += "%-15.3f "%v[1]
      else: 
        if (n < minEventNumber) and (not ('lumi' in s)) and (not ('scale' in s)): 
          l += "%-15s "%"-"
        else:
          vstr = "%.3f/%.3f"%(v[0],v[1])
          if s == "CMS_hgg_SigmaEOverEShift":
            print(vstr, p, c, n)
          l += "%-15s "%vstr
    return l
  else:
    print " --> [ERROR] systematic %s: value does not have type string or list for (%s,%s). Leaving..."%(s['title'],p,c)
    sys.exit(1)

def writeMCStatUncertainty(f,d,options):
 
  # Remove all rows from dataFrame with prune=1 (includes NoTag)
  if options.prune:
    mask = (d['prune']==0)
    d = d[mask]

  # Separate nuisance for each cat * years: ~Barlow-Beeston-Lite approach
  for scat in d.cat.unique():
    for year in options.years.split(","):
      # Extract sval
      mask = (d['year']==year)&(d['cat']==scat)&(d['type']=='sig')
      sumw = d[mask]['nominal_yield'].sum()
      sumw2 = d[mask]['sumw2'].sum()
      scval = [1+(math.sqrt(sumw2)/sumw)]
      d[d['type']=='sig']
      stitle = "MCStat_%s_%s"%(year,scat)
      sprior = "lnN"
      lsyst = '%-50s  %-10s    '%(stitle,sprior)
      # Loop over categories and then iterate over rows in category
      for cat in d.cat.unique():
        for ir,r in d[d['cat']==cat].iterrows():
          if r['proc'] == "data_obs": continue
          if (not (year in r['proc'])) or (r['proc'] == "bkg_mass"): 
            # Do not assign the MC uncertainty to procs not matching the year/era and to the data
            lsyst += "%-15s "%"-"
            continue
          sval = scval if cat == scat else '-'
          # Extract value and add to line (with checks)
          lsyst = addSyst(lsyst,sval,stitle,r['proc'],cat,r['numEvents'])
      # Remove final space from line and add to file
      f.write("%s\n"%lsyst[:-1])
  
  return True

def writePdfIndex(f,d,options):
  f.write("\n")
  years = list(set([re.search(r'\d{4}', item).group() for item in options.years.split(",") if re.search(r'\d{4}', item)]))
  for cat in d[~d['cat'].str.contains("NOTAG")].cat.unique():
    if options.mergeYears:
      indexStr = "pdfindex_%s_13TeV"%cat
      f.write("%-55s  discrete\n"%indexStr)
    else:
      for year in years:
        indexStr = "pdfindex_%s_%s_13TeV"%(cat, year)
        f.write("%-55s  discrete\n"%indexStr)
  return True

def writeBreak(f):
  lbreak = '----------------------------------------------------------------------------------------------------------------------------------'
  f.write("%s\n"%lbreak)

