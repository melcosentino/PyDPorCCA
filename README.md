 #D-PorCCA
  D-PorCCA is a standalone desktop application aiming at providing researchers with summary data to monitor harbour
   porpoises from continuous recordings.  
   D = Detector
   PorCC = Porpoise Click Classifier
   A = Application

#### Detector
The Detector is an adaptation of PAMGuard's Click Detector Module (translated by Clea Parcerisas) that uses a trigger filter selecting clips that can 
potentially be porpoise clicks (i.e., high energy content between 100 and 150 kHz)

#### PorCC: 
> Cosentino, M., Guarato, F., Tougaard, J., Nairn, D., Jackson, J. C., & Windmill, J. F. C. (2019). 
> Porpoise click classifier (PorCC): A high-accuracy classifier to study harbour porpoises (*Phocoena phocoena*) in the wild . 
> The Journal of the Acoustical Society of America, 145(6), 3427–3434. https://doi.org/10.1121/1.5110908

In DPorCCA, PorCC is implemented via de package pyporcc, developed by Clea Parcerisas

## Note
DPorCCA is still under development


 ## INFORMATION ABOUT PYTHON FILES
- click_trains.py: a series of functions to identify and classify click trains
- isoutlier.py: identifies outliers in the data - translated from Matlab by Clea Parcerisas
- sunrise.py: estimates sunset and sunrise times for a given location in a given day (GMT) from latitude and longitude
- create_settings_file.py: generates setting files when using the Detector and PorCC.

### Examples:


## Outputs
CP = .csv file 
Fields: 
- *id*: identification number
- *date*: date as dd-mmm-yyyy hh:mm:ss 
- *start_sample*: sample where the signal begins  
- *duration*: duration estimated as the 80% of the energy of the signal 
- *CF*: centroid frequency 
- *BW*: the -3dB bandwidth
- *ratio*: ratio between the peak and centroid frequency
- *XC*: maximum value of the cross-correlation coefficient carried out against a typical porpoise click, 
- *Q*: defined as the RMS bandwidth divided the centroid frequency
- *pyporcc*: class assigned by PorCC (1: high-quality click, 2: low-quality click)


Clicks = .csv file 
Fields: 
- *id*: identification number
- *date*: date as dd-mmm-yyyy hh:mm:ss 
- *start_sample*: sample where the signal begins  
- *duration*: duration estimated as the 80% of the energy of the signal 
- *CF*: centroid frequency 
- *BW*: the -3dB bandwidth
- *ratio*: ratio between the peak and centroid frequency
- *XC*: maximum value of the cross-correlation coefficient carried out against a typical porpoise click, 
- *Q*: defined as the RMS bandwidth divided the centroid frequency
- *pyporcc*: class assigned by PorCC (1: high-quality click, 2: low-quality click)
- *ICI*: inter-click interval (ms)
- *CPS*: clicks per second

CTInfo = .csv file 
Fields: 
- *CTNum*: click train number
- *Date*: date as dd-mmm-yyyy hh:mm:ss
- *DayNight*: describes whether the event occurred during the day or at night  
- *Length*: number of clicks  
- *CTType*: click train type 
- *Behav*: behaviour of the animal  
- *Calf*: 
- *Notes*: 






























