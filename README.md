 #D-PorCCA
  D-PorCCA is a standalone desktop application with the goal of providing researchers with summary data for the 
   monitoring of harbour porpoises.  
   D = Detector
   PorCC = Porpoise Click Classifier
   A = Application

#### Detector
D-PorCCA includes a detector for impulsive sounds, regularly called "click detector" and it is an adaptation of 
PAMGuard's Click Detector Module (translated by Clea Parcerisas), which uses a trigger filter to select potential
porpoise clicks (i.e., high energy content between 100 and 150 kHz)

#### PorCC: 
The porpoise click classifier (PorCC) in D-PorCCA was developed by Cosentino et al and separates clicks into either of 
three categories: high-frequency noise, high-quality porpoise click, low-quality porpoise click. More information can be 
found here:
> Cosentino, M., Guarato, F., Tougaard, J., Nairn, D., Jackson, J. C., & Windmill, J. F. C. (2019). 
> Porpoise click classifier (PorCC): A high-accuracy classifier to study harbour porpoises (*Phocoena phocoena*) in the wild . 
> The Journal of the Acoustical Society of America, 145(6), 3427–3434. https://doi.org/10.1121/1.5110908

PorCC is implemented via the package pyporcc, developed by Clea Parcerisas.

## Note
DPorCCA is still under development


 ## FILES INFORMATION 
- click_trains.py: functions to identify, clean and classify click trains
- isoutlier.py: identifies outliers in the data - function translated from Matlab by Clea Parcerisas
- sunrise.py: estimates sunset and sunrise times for a given location on a given day (GMT) from latitude and longitude
- create_settings_file.py: generates setting files when using the Detector and PorCC.
- GUIDPorCCA.py: main script where the app is created

- requirements.txt: list of packages needed for the proper functioning of D-PorCCA. During installation, the packages 
listed in this file are installed automatically.

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
- *Calf*: *NOT IMPLEMENTED YET*
- *Notes*: 






























