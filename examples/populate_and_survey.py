#!/usr/bin/python

# import the populate and dosurvey modules
from populate import Populate
from dosurvey import DoSurvey

# run the Populate.generate script, setting parameters as necessary.
# the resulting model is stored in the variable 'pop.pop'
# Any unspecified variables use the default values (see populate.py
# for more details)

pop = Populate()
pop.generate(1038, 
             surveyList=['PMSURV'],
             radialDistType='lfl06',
             siDistPars=[-1.41, 0.96], # non-standard SI distribution
             duty=6.,
             electronModel='lmt85',
             nostdout=True # switches off output to stdout
             )

# now run "dosurvey" on the model. Provide a list of surveys to use
ds = DoSurvey(popmodel = pop.pop)
ds.run(['PMSURV'], nostdout=True)

# write out the survey results however you like
ds.write(nores=False, summary=False, asc=False)

# write out the population model, if required
pop.write(outf="populate.model")
