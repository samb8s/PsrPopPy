#!/bin/tcsh

foreach i (`seq 150`)
    dosurvey -f populate.model -surveys LOFAR | grep detected >> noscint.list
end

# 1151 +- 17

foreach i (`seq 150`)
    dosurvey -f populate.model -surveys LOFAR --scint | grep detected >> scint.list
end

# 
