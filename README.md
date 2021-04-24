 #D-PorCCA
  D-PorCCA is a standalone desktop application with the goal to provide researchers with summary data to monitor harbour
   porpoises from continuous recordings.  
   D = Detector
   PorCC = Porpoise Click Classifier 
   A = Application
  
Everything about PorCC can be found here: 
> Cosentino, M., Guarato, F., Tougaard, J., Nairn, D., Jackson, J. C., & Windmill, J. F. C. (2019). 
> Porpoise click classifier (PorCC): A high-accuracy classifier to study harbour porpoises ( Phocoena phocoena ) in the wild . 
> The Journal of the Acoustical Society of America, 145(6), 3427–3434. https://doi.org/10.1121/1.5110908

In DPorCCA it is implemented via de package PyPorCC developed by Clea Parcerisas

The PyPorCC package also provides an adapted alternative to PAMGuard's click detector that using a filter and a trigger 
function selects clips that can potentially be clicks (high enery in the right frequency band)

## Note
This application is still under development


 ## INFORMATION ABOUT PYTHON FILES
- click_trains.py: a series of functions to identify and classify click trains
- isoutlier.py: identifies outliers in the data - from Matlab
- sunrise.py: estimates sunset and sunrise times for a given location in a given day
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
- *ICI*: inter-click interval
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






























