#!/usr/bin/python

# import the populate and dosurvey modules
import populate
import dosurvey

# run the Populate.generate script, setting parameters as necessary.
# the resulting model is stored in the variable 'pop'
# Any unspecified variables use the default values (see populate
# for more details)

pop = populate.generate(1038, 
               surveyList=['PMSURV'],
               radialDistType='lfl06',
               siDistPars=[-1.41, 0.96], # non-standard SI distribution
               duty_percent=6.,
               electronModel='lmt85',
               nostdout=True # switches off output to stdout
               )

# now run "dosurvey.run" on the model. Provide a list of surveys to use
surveyPopulations = dosurvey.run(pop, ['PMSURV'], nostdout=True)

# write out the survey results however you like
dosurvey.write(surveyPopulations, nores=False, summary=False, asc=False)

# write out the population model, if required
pop.write(outf="populate.model")
