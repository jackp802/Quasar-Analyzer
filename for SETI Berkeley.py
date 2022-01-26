# -*- coding: utf-8 -*-

from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import make_interp_spline
from astropy.modeling import models, fitting
from astropy import units as u
import os
from matplotlib.backends.backend_pdf import PdfPages
from scipy.integrate import quad

cwd = os.getcwd()
directory=r"C:\Users\Owner\Desktop\Fits\BAL"
graphs=PdfPages('GRAPHS.pdf')
normalgraphs=PdfPages('NGRAPHS.pdf')

c4fit=np.array([])
masterpowerlaw=np.array([])
gaussianfit=np.array([])
c4lflux=np.array([])

totalfiles=0
for filename in os.listdir(directory):
    if filename.endswith(".fits"):
        print("Name of file: ", filename)
        print('File number: ',totalfiles)
        from astropy.io import fits
                         
                #=os.path.join(r'C:\Users\Owner\Desktop\Fits\BAL',filename)
                #print(nam)
                #name=fits.util.get_testdata_filepath(nam)
        name = os.path.join(directory, filename)
        hdul=fits.open(name)
        # when there is # in the front of code, just means that I am not
        #currently using that code to make the output prettier and more simple
        plt.figure()
        totalfiles+=1
        with fits.open(name) as hdul:
                   # hdul.info()
                #print("====================HEADER=====================================")
                y=hdul[0].header
                #print(y)
                hdul.close()
                hdul=fits.open(name)
                #print("========================hdul[x].columns (maybe header too?)========")
                
                hdu1=hdul[1]
                #print(hdul[1].columns)
                cols=hdul[1].columns
                f=cols['flux']
                #print("====================Wavelength,flux,redshift=======================")
                
                
                data=hdul[1].data
                data2=hdul[2].data
                z=data2['z']
                flux=data['flux']
                
                
                logwave=data['loglam']
                Rwave=10**logwave
                #print(Rwave)
                #print(flux)
                #print(z)
                #print("==============================================================")
                #plt.plot(Rwave,flux)
                
                wave=Rwave/(1+z)
                #plt.show()
                #print("Above is not corrected for redshift")
                #print(' ')
                
                                                    #1+Z = loglam_obs/loglam_emit
                
                c=3*(10**5)
                dlam=(1-(1/(1+z)))*(1+z)
                recvel=c*dlam
                hubb=73
                dmpc=recvel/hubb
                d=dmpc*10**6
                dm=d*3.086*10**18
                #print('Recessive Velocity', recvel, 'km/s')
                #print('Distance',d,'pc')
                import math
                
                
                arr=np.where(np.logical_and(wave>=1900,wave<=2350))
                
                cf=sum(flux[arr]/len(flux[arr]))
                stdev=np.std(flux[arr])
                #print('STDEV:',stdev)
                #cf=np.interp(2200, wave,flux)
                cferror=stdev/(len(flux[arr])-1)
                #print('Continuum Flux:',cf)
                #print('Its error:',cferror*100)
                
                L=cf*4*math.pi*(dm**2)
                #print('Luminosity',L,'erg/s')
                
                
                #===========================Plots====================
                #==================Smoothed/normal
                #print('Flux vs Emitted Wavelength')
                xnew=np.linspace(wave.min(),wave.max(),500)
                
                plt.title("Flux vs Wave with Fitted Red Line")
                plt.xlim()
                plt.ylim(-1,15)
                plt.xlabel("Wavelength")
                plt.ylabel("Flux")
                spline=make_interp_spline(wave,flux)
                ynew=spline(xnew)
                #plt.plot(xnew,ynew)
                #=====================Fitted line
                from scipy.optimize import curve_fit
                from numpy import arange
                def powerlaw(x,a,b):
                    return a*x**b
                
                fit1,_=curve_fit(powerlaw,wave,flux)
                a,b=fit1
                xfit=arange(wave[0],wave[-1])
                yfit=powerlaw(xfit,a,b)
                #plt.plot(xfit,yfit,color='red')
                #plt.show()
                #===================Chi^2
                nflux=np.array([])
                for p in wave:
                    nflux=np.append(nflux,a*p**b)
                from scipy.stats import chisquare
                #print(chisquare(flux,nflux))
                #print('===================')
                #===================Fitted line with continuum subtracted
                gflux=np.subtract(flux,nflux)
                spline=make_interp_spline(wave, gflux)
                ynew=spline(xnew)
                #plt.plot(xnew,ynew)
                #plt.plot(wave,flux*0,color='red', label='X-Axis')
                #plt.title("Flux vs Wave with Continuum Flux subtracted")
                #plt.xlabel("Wavelength")
                #plt.ylabel("Flux")
                plt.ylim(-1,10)
                #plt.legend()
                #plt.show()
                
                #Below is locoation to save as png
                #plt.savefig('C:/Users/Owner/desktop/Spectra/waveflux.png')
                
                #================Plot with Limit (for zooming in)
                #print('Plot with limits')
                spline=make_interp_spline(wave,flux)
                ynew=spline(xnew)
                
                arr=np.where(np.logical_and(xnew>=1200,xnew<=2100))
                newwave=xnew[arr]
                newflux=ynew[arr]
                #plt.plot(newwave,newflux)
                
                #plt.plot(xnew,ynew*0)
                plt.xlim(1200,1800) #Carbon 3. Fit the contiuum. subtract from bump. integrate bump
                plt.ylim(-1,6)
                #=================Emission Line Fitting
                arr=np.where(np.logical_and(xnew>=1300,xnew<=1400))
                newwave=xnew[arr]
                newflux=ynew[arr]
                arr=np.where(np.logical_and(xnew>=1600,xnew<=1800))
                x=xnew[arr]
                y=ynew[arr]
                x1=np.concatenate((newwave,x))
                y1=np.concatenate((newflux,y))
                
                fit=fitting.LinearLSQFitter() 
                line_init=models.Polynomial1D(degree=3)
                fitted=fit(line_init,x1,y1)
                #plt.plot(x1,fitted(x1),label='fit')
                
                #plt.legend()
                #plt.show()
                #===================================
                #=======================Fit subtracted from Emission Line
                from specutils.spectra import Spectrum1D
                from specutils.fitting import fit_lines
                arr=np.where(np.logical_and(xnew>=1300,xnew<=1800))
                newflux=ynew[arr]
                newwave=xnew[arr]
                gflux=np.subtract(newflux,fitted(newwave))
                #plt.plot(newwave,gflux)
                plt.xlim(1300,1700)
                plt.ylim(-1,2)
                #plt.show()
                #==============================================Gaussian Attempt
                arr=np.where(np.logical_and(xnew>=1350,xnew<=1425))
                newwave=xnew[arr]
                newflux=ynew[arr]
                arr=np.where(np.logical_and(xnew>=1525,xnew<=1650))
                x=xnew[arr]
                y=ynew[arr]
                
                newwave=np.concatenate((newwave,x))
                newflux=np.concatenate((newflux,y))
                gflux=np.subtract(newflux,fitted(newwave))
                spectrum=Spectrum1D(flux=gflux*u.Jy,spectral_axis=newwave*u.AA)
                g_init = models.Gaussian1D(amplitude=2.5*u.Jy, mean=1549*u.AA, stddev=40.*u.AA)
                g_init.stddev.fixed=False
                g_init.amplitude.fixed=False
                g_init.mean.fixed=False
                
                
                
                fit_g = fit_lines(spectrum,g_init)
                arr=np.where(np.logical_and(xnew>=1300,xnew<=1650))
                newwave=xnew[arr]
                newflux=ynew[arr]
                gflux=np.subtract(newflux,fitted(newwave))
                #plt.plot(newwave,fit_g(newwave*u.AA),label='GaussFit',color='red')
                #plt.plot(xnew,ynew*0)
                plt.xlim(1300,1700)
                plt.ylim(-1,3)
                #plt.show()
                
                
                
                
                arr=np.where(np.logical_and(xnew>=1450,xnew<=1700))
                newwave=xnew[arr]
                newflux=ynew[arr]
                gflux=np.subtract(newflux,fitted(newwave))
                max1=np.max(gflux)
                maxi=np.where((gflux==max1))
                wav=newwave[maxi]
                
                arr2=np.where(np.logical_and(xnew>=wav,xnew<=1700))
                farr2=np.flip(arr2)
                barr2=np.concatenate((farr2[0],arr2[0]))
                larr=int(len(barr2))
                
                
                arr=np.where(np.logical_and(xnew>=1300,xnew<=1700))
                
                arr=arr[0][1:int(len(arr))-larr]
                arrf=np.concatenate((arr,barr2))
                arr=np.where(np.logical_and(xnew>=1300,xnew<=1700))
                newwave=xnew[arr]
                newflux=ynew[arrf]
                
                marr=np.where(np.logical_and(xnew>=1300,xnew<=1700))
                newwave=xnew[arr]
                newflux=ynew[arrf]
                
                mflux=np.subtract(newflux,fitted(newwave))
                
                plt.plot(xnew,ynew*0)
                mwave=xnew[marr]
                plt.plot(newwave,mflux,color='blue')
               
                #Gaussian 
                gflux=np.subtract(newflux,fitted(newwave))
                spectrum=Spectrum1D(flux=gflux*u.Jy,spectral_axis=newwave*u.AA)
                g_init = models.Gaussian1D(amplitude=np.amax(gflux)*u.Jy, mean=1549*u.AA, stddev=40.*u.AA)
                g_init.stddev.fixed=False
                g_init.amplitude.fixed=False
                g_init.mean.fixed=True               
                
                fit_g = fit_lines(spectrum,g_init)

                arr=np.where(np.logical_and(xnew>=1350,xnew<=1650))
                newwave=xnew[arr]
                newflux=ynew[arr]
                gflux=np.subtract(newflux,fitted(newwave))
                plt.plot(newwave,fit_g(newwave*u.AA),label='GaussFit',color='red')
                plt.plot(xnew,ynew*0,color='orange')
                plt.xlim(1300,1700)
                plt.ylim(-1,10)
                arr=np.where(np.logical_and(xnew>=1300,xnew<=1800))
                newflux=ynew[arr]
                newwave=xnew[arr]
                gflux=np.subtract(newflux,fitted(newwave))
                plt.plot(newwave,gflux, color='green')
                
                
                #Below saves the graphs as a PDF
                #plt.savefig(graphs,format='pdf')
                plt.title("Carbon4 Absoprtion line Fit")
                plt.show()
                print(""" Green line is actual data. Blue line is data that is mirrored
from the center of C4's emission line (1549nm). Red line is the gaussian
that is fitted over it.
                      """)
                #======================End of fitted curve
                #Start of normalized data
                arr=np.where(np.logical_and(xnew>=1450,xnew<=1600))
                newflux=ynew[arr]
                newwave=xnew[arr]
                #Gauss fit from 1450-1600
                gaussf=np.array([fit_g(newwave*u.AA)]) #Converting gaussian function into array
                
                
                
                arr=np.where(np.logical_and(xnew>=1200,xnew<=1450))
                x2=xnew[arr]
                y2=ynew[arr]
                arr=np.where(np.logical_and(xnew>=1600,xnew<=1800))
                x3=xnew[arr]
                y3=ynew[arr]
                
                
                #Creating fit from 1200-1450
                fit2=fitting.LinearLSQFitter() 
                line_init2=models.Polynomial1D(degree=3)
                fitted2=fit(line_init2,x2,y2)
                fitmean2=np.mean(fitted2(x2))
                
                #Creating fit from 1600-1800
                fit3=fitting.LinearLSQFitter() 
                line_init3=models.Polynomial1D(degree=3)
                fitted3=fit(line_init2,x3,y3)
                fitmean3=np.mean(fitted3(x3))
                
                #Getting avg value from the fitted to add to gaussian so it is on same y-value approx
                fitavg=(fitmean2+fitmean3)/2
                
                normg=gaussf + fitavg  #Combining all the arrays THESE ARE THE VALUES FOR FLUX
                norm1=np.concatenate((fitted2(x2),normg[0]))
                norm2=np.concatenate((norm1,fitted3(x3)))
                #Norm2 is the fitted line with the gaussian ontop of the emission line
                
                arr=np.where(np.logical_and(xnew>=1200,xnew<=1800))
                newwave=xnew[arr]
                newflux=ynew[arr]
                NORMALIZED=np.divide(newflux,norm2)
                plt.plot(newwave,NORMALIZED)
                plt.title('NORMALIZED')
                plt.savefig(normalgraphs,format='pdf')
                plt.plot(newwave,np.ones(len(newwave)),color='black')
                plt.ylim(0,2)
                plt.show()
                #############################Print value of min flux's wavelength
                
                mflux  #[arrf] Mirror flux non guss
                mwave  #[marr] Normal wavelength 1300-1700
                #fit_g(newwave*u.AA) Gauss between 1350-1650
                
                #Below is finding the lowest flux point
                
                arr=np.where(np.logical_and(xnew>=1400,xnew<=1600))
                newwave=xnew[arr]
                newflux=ynew[arr]
                gflux=np.subtract(newflux,fitted(newwave))
                Lflux=np.where(gflux==np.amin(gflux))
                print("Lowest Flux Wavelength: ",newflux[Lflux])
                
                c4lflux=np.append(c4lflux,Lflux)
                masterpowerlaw=np.append(masterpowerlaw,[yfit])
                c4fit=np.append(c4fit,[fitted])
                gaussianfit=np.append(gaussianfit,[fit_g])
                
                #Wavelength to velocity
                velo=np.array([])
                for x in xnew:
                    vel=((x-1549)/1549)*3*10**5
                    velo=np.append(velo,vel)
                    
                
                #Below shows the plot of flux vs velocity
                #plt.xlim(-50000,50000)
                #plt.title('velo')
                #plt.plot(velo,ynew)
                #plt.show()
                arr=np.where(np.logical_and(xnew>=1200,xnew<=1800))
                newwave=xnew[arr]
                newflux=ynew[arr]
                
                
                cv=np.array([])
                length=3
                c=1
                #1200-1800
                
                    #np.ones(200) np.zeros(200)
                
                
                velo1=np.array([])
                for x in xnew[arr]:
                    vel=((x-1549)/1549)*3*10**5
                    velo1=np.append(velo1,vel)
                    
                velarr=np.where(np.logical_and(velo1>=-25000,velo1<=3000))
                
                
                #Below is currently being used to find the Balnicity index
                c=1
                for x in NORMALIZED[velarr]:
                    if x<0.9 and c<length:
                        c+=1
                        
                    elif x>0.9 and c>=length:
                        
                        cv=np.append(cv,np.ones(c))
                        c=1
                    elif x>0.9:
                        
                        cv=np.append(cv,np.zeros(c))
                        c=1
                
                
                
                def BI(x,a,b):
                    return -(1-(a/0.9))*b
                
                a=NORMALIZED[velarr]
                b=cv
                #BIinput=BI(a,b)
                BIarr=np.array([])
                BIarray=np.array([])
                for a in NORMALIZED[velarr]:
                    for b in cv:
                        
                        BIarr=np.append(BIarr,a*2.1*b)
                        #above two make it for every cv and flux value, integrate this.
                        ans=quad(BI,-25000,3000,args=(a,b))
                        #above is the integration code
                        #below is appending to an array
                        BIarray=np.append(BIarray,ans)
                            
                #print(np.sum(BIarr))
                #print(BIarray)
                #print(len(cv),'vs',len(NORMALIZED[velarr]))            
                
                
                
                #This is plotted the smoothed spectra
                plt.plot(xnew,ynew,color='r')
                plt.title('Flux vs Wavelenth, smoothed')
                plt.xlabel('Wave')
                plt.ylabel('Flux')
                plt.ylim(0,10)
                plt.xlim(1250,3000)
                plt.show()
                #Use below to show raw data plot
                #plt.ylim(0,10)
                #plt.xlim(1250,3000)
                #plt.plot(wave,flux)
                #plt.show()
                
                
                
                
                print('================================================')
                hdul.close()

print(totalfiles, "THIS IS THE TOTAL AMOUNT OF FILES RAN")
graphs.close()
normalgraphs.close()

