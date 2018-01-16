#!/bin/bash

# Script for running several VisIt python scripts

# Check to see if all files are present first!! 

/Applications/VisIt.app/Contents/Resources/bin/visit -nowin -cli -s aorta_P.py
mv *.curve curvefiles/
/Applications/VisIt.app/Contents/Resources/bin/visit -nowin -cli -s aorta_Um.py
mv *.curve curvefiles/
/Applications/VisIt.app/Contents/Resources/bin/visit -nowin -cli -s aorta_Ux.py
mv *.curve curvefiles/
/Applications/VisIt.app/Contents/Resources/bin/visit -nowin -cli -s connect_Um.py
mv *.curve curvefiles/
/Applications/VisIt.app/Contents/Resources/bin/visit -nowin -cli -s connect_Uy.py
mv *.curve curvefiles/
/Applications/VisIt.app/Contents/Resources/bin/visit -nowin -cli -s topcross_Um.py
mv *.curve curvefiles/
/Applications/VisIt.app/Contents/Resources/bin/visit -nowin -cli -s topcross_Ux.py
mv *.curve curvefiles/
/Applications/VisIt.app/Contents/Resources/bin/visit -nowin -cli -s vena_P.py
mv *.curve curvefiles/
/Applications/VisIt.app/Contents/Resources/bin/visit -nowin -cli -s vena_Um.py
mv *.curve curvefiles/
/Applications/VisIt.app/Contents/Resources/bin/visit -nowin -cli -s vena_Ux.py
mv *.curve curvefiles/


/Applications/VisIt.app/Contents/Resources/bin/visit -nowin -cli -s profile_topcross_Um.py