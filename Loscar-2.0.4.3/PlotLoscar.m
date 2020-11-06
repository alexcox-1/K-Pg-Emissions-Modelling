%================================================================%
%
% PlotLoscar.m (Matlab script)
%
% LOSCAR Model: Long-term Ocean-atmosphere-Sediment 
%                 CArbon cycle Reservoir Model
% 
% *** LOSCAR comes with ABSOLUTELY NO WARRANTY ***		 
% *** Use at your own risk. DO NOT DISTRIBUTE  ***
%
% When using LOSCAR, cite as:
%
% Zeebe, R. E., Zachos, J. C., and Dickens, G. R. Carbon dioxide
% forcing alone insufficient to explain Paleocene-Eocene Thermal
% Maximum warming. Nature Geoscience, doi:10.1038/NGEO578 (2009) 
% 
% Richard E. Zeebe
% School of Ocean and Earth 
% Science and Technology 
% Department of Oceanography 
% University of Hawaii at Manoa
% 1000 Pope Road, MSB 504
% Honolulu, HI 96822, USA
% email: loscar.model@gmail.com
%
%
% Purpose: Matlab script to plot Loscar model output.
%
% Matlab Version: R2007b and up.
%
% updates:
%
%  09/03/18 removed location '4' from legends
%  10/17/11 include Tethys
%  09/16/11 include d13C
%  09/11/11 include dissolved O2 (dox).
%  04/06/11 new file
%
%================================================================%

clear all all

% 0: load data only, 1: load and plot
plotflag = 1;

% defaults: don't plot unless .dat files exist
fchem  = 0;
ftemp  = 0;  
fdox   = 0;
f13c   = 0;
fsed   = 0; 
f13csd = 0; 
ftys   = 0; 

% plot CO2 chemistry ?
if(exist('co3.dat','file') ~= 0)
fchem = 1;
end
% plot temperature ?
if(exist('tcb.dat','file') ~= 0)
ftemp = 1;
end
% plot dissolved O2 ?
if(exist('dox.dat','file') ~= 0)
fdox  = 1;
end
% plot 13C ?
if(exist('dicc.dat','file') ~= 0)
f13c  = 1;
end
% plot sediment variables ?
if(exist('fca.dat','file') ~= 0)
fsed  = 1; 
end
% plot 13C sediment variables ?
%if(exist('fcca.dat','file') ~= 0)
%f13csd  = 1; 
%end
% plot Tethys ?
if(exist('ccdt.dat','file') ~= 0)
ftys  = 1; 
end


km2m = 1.e3;

if(ftys)
Nb = 13;
else
Nb = 10;
end
Ns = 13;


%============== load model output ====================%

% set directory ('' = current)
dirstr = '';

tv    = load([dirstr      'tmv.dat']);
dic   = load([dirstr      'dic.dat']);
alk   = load([dirstr      'alk.dat']);
po4   = load([dirstr      'po4.dat']);
pco2a = load([dirstr    'pco2a.dat']);
if(fchem)
co3   = load([dirstr      'co3.dat']);
ph    = load([dirstr       'ph.dat']);
omclc = load([dirstr 'omegaclc.dat']);
omarg = load([dirstr 'omegaarg.dat']);
end
if(ftemp)
tcb   = load([dirstr      'tcb.dat']);
end
if(fdox)
dox   = load([dirstr      'dox.dat']);
end
if(f13c)
dicc  = load([dirstr     'dicc.dat']);  
d13c  = load([dirstr     'd13c.dat']);
d13ca = load([dirstr    'd13ca.dat']);
end
if(fsed)
dsv   = load([dirstr      'dsv.out']);
fca   = load([dirstr      'fca.dat']);
fci   = load([dirstr      'fci.dat']);
fcp   = load([dirstr      'fcp.dat']);
ccda  = load([dirstr     'ccda.dat'])/km2m;
ccdi  = load([dirstr     'ccdi.dat'])/km2m;
ccdp  = load([dirstr     'ccdp.dat'])/km2m;
end
if(f13csd)
fcca  = load([dirstr     'fcca.dat']);
fcci  = load([dirstr     'fcci.dat']);
fccp  = load([dirstr     'fccp.dat']);
end
if(ftys)
fct   = load([dirstr      'fct.dat']);    
ccdt  = load([dirstr     'ccdt.dat'])/km2m;
  if(f13csd) 
    fcct  = load([dirstr     'fcct.dat']);
  end
end


lt = length(tv);

% print pCO2 values
pCO2_max_pCO2_final = [max(pco2a) pco2a(lt)];
%format long
pCO2_max_pCO2_final
%format short



if(plotflag)
%=============== plot parameters =====================%
i = 1;
fs  = 14;
lfs = 10;
if(ftys)
kkv   = [1 2 3 10 11];
else
kkv   = [1 2 3 10];
end
xlm   = [tv(1) tv(lt)];
ylmfc = [-10. 100.];
cs    = 'gggkkkrrrbgkr';
css   = 'bgrmckybgrmck';
sstr  = '- ---.- ---.- ---. -: : : ';
lstr0 = 'LALILPIAIIIPDADIDP HLTITDT';
for k=1:Nb
lstr(k,:) = sprintf('%s',lstr0(2*k-1:2*k));
end
for k=1:length(kkv)
lstrs(k,:) = sprintf('%s',lstr0(2*kkv(k)-1:2*kkv(k)));
end

if(fsed)
for l=1:Ns
dstr(l,:) = sprintf('%4.0f m',dsv(l));
end;
if(ftys)
ccdstr = ['Atl ';'Ind ';'Pac ';'Teth'];
else
ccdstr = ['Atl';'Ind';'Pac'];
end
ocncol = 'bmkc';
end

%===================== plot ==========================%

%--------------------- dic ---------------------------%
figure(i); i=i+1;
clf 
box  on
hold on
for k=1:Nb
plot(tv,dic(:,k),sstr(2*k-1:2*k),'Color',cs(k))
end
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
xlabel('Time (y)')
ylabel('DIC (mmol kg^{-1})')
Hl=legend(lstr);
set(Hl,'FontSize',lfs)

%---------------------- alk --------------------------%
figure(i); i=i+1;
clf 
box  on
hold on
for k=1:Nb
plot(tv,alk(:,k),sstr(2*k-1:2*k),'Color',cs(k))
end
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
xlabel('Time (y)')
ylabel('TA (mmol kg^{-1})')
Hl=legend(lstr);
set(Hl,'FontSize',lfs)

%---------------------- po4 --------------------------%
figure(i); i=i+1;
clf 
box  on
hold on
for k=1:Nb
plot(tv,po4(:,k),sstr(2*k-1:2*k),'Color',cs(k))
end
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
xlabel('Time (y)')
ylabel('PO_4 (\mumol kg^{-1})')
Hl=legend(lstr);
set(Hl,'FontSize',lfs)

if(ftemp)
%----------------- temperature -----------------------%
figure(i); i=i+1;
clf 
box  on
hold on
for k=1:Nb
plot(tv,tcb(:,k),sstr(2*k-1:2*k),'Color',cs(k))
end
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
xlabel('Time (y)')
ylabel('Temperature (deg C)')
Hl=legend(lstr);
set(Hl,'FontSize',lfs)

end % ftemp

if(fdox)
%----------------- dissolved O2 ----------------------%
figure(i); i=i+1;
clf 
box  on
hold on
for k=1:Nb
plot(tv,dox(:,k),sstr(2*k-1:2*k),'Color',cs(k))
end
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
xlabel('Time (y)')
ylabel('Dissolved O_2 (mol m^{-3})')
Hl=legend(lstr);
set(Hl,'FontSize',lfs)

end % ftemp

%--------------------- 13C -------------------------%
figure(i); i=i+1;
clf 
box  on
hold on
for k=1:Nb
plot(tv,d13c(:,k),sstr(2*k-1:2*k),'Color',cs(k))
end
Dam   = 7.0;
plot(tv,d13ca+Dam,'m')
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
set(gca,'YDir','r')
xlabel('Time (y)')
ylabel('\delta^{13}C (^o/oo)')
Hl=legend([lstr;'At']);
set(Hl,'FontSize',lfs)
str = sprintf('At = Atm + %.1f^o/oo',Dam);
Ht=text(.5,.1,str,'Units','n');
set(Ht,'FontSize',lfs,'Hor','c')

%---------------------- pCO2 -------------------------%
figure(i); i=i+1;
plot(tv,pco2a,'b-')
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
xlabel('Time (y)')
ylabel('Atmospheric CO_2 (ppmv)')

if(fchem)
%---------------------- CO3 --------------------------%
figure(i); i=i+1;
clf 
box  on
hold on
for k=1:Nb
plot(tv,co3(:,k),sstr(2*k-1:2*k),'Color',cs(k))
end
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
xlabel('Time (y)')
ylabel('[CO_3^{2-}] (\mumol kg^{-1})')
Hl=legend(lstr);
set(Hl,'FontSize',lfs)

%----------------------  pH --------------------------%
figure(i); i=i+1;
clf 
box  on
hold on
for k=1:Nb
plot(tv,ph(:,k),sstr(2*k-1:2*k),'Color',cs(k))
end
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
xlabel('Time (y)')
ylabel('pH');
Hl=legend(lstr);
set(Hl,'FontSize',lfs)

%---------------- Omega-calcite ----------------------%
figure(i); i=i+1;
clf 
box  on
hold on
for k=kkv
plot(tv,omclc(:,k),sstr(2*k-1:2*k),'Color',cs(k))
end
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
xlabel('Time (y)')
ylabel('\Omega-calcite')
Hl=legend(lstrs);
set(Hl,'FontSize',lfs)

%---------------- Omega-aragonite --------------------%
figure(i); i=i+1;
clf 
box  on
hold on
for k=kkv
plot(tv,omarg(:,k),sstr(2*k-1:2*k),'Color',cs(k))
end
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
xlabel('Time (y)')
ylabel('\Omega-aragonite')
Hl=legend(lstrs);
set(Hl,'FontSize',lfs)

end % fchem

if(fsed)
%-------------------- fc vs. depth -------------------%
figure(i); i=i+1;
clf 
box  on
hold on    
plot(fca(lt,:)*1e2,dsv/km2m,'-d','Color',ocncol(1));
plot(fci(lt,:)*1e2,dsv/km2m,'-s','Color',ocncol(2));
plot(fcp(lt,:)*1e2,dsv/km2m,'-o','Color',ocncol(3));
if(ftys)
plot(fct(lt,:)*1e2,dsv/km2m,'-o','Color',ocncol(4));
end
hold off
set(gca,'FontSize',fs)
set(gca,'YDir','reverse')
xlabel('CaCO_3 (wt %)')
ylabel('Depth (km)')
Hl=legend(ccdstr);
set(Hl,'FontSize',lfs)

%----------------------- fca -------------------------%
figure(i); i=i+1;
clf 
box  on
hold on
for l=1:Ns
plot(tv,fca(:,l)*1e2,'--','Color',css(l))
end
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
set(gca,'YLim',ylmfc)
title('Atlantic Sediments')
xlabel('Time (y)')
ylabel('CaCO_3 (wt %)')
Hl=legend(dstr);
set(Hl,'FontSize',lfs)

%----------------------- fci -------------------------%
figure(i); i=i+1;
clf 
box  on
hold on
for l=1:Ns
plot(tv,fci(:,l)*1e2,'--','Color',css(l))
end
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
set(gca,'YLim',ylmfc)
title('Indian Sediments')
xlabel('Time (y)')
ylabel('CaCO_3 (wt %)')
Hl=legend(dstr);
set(Hl,'FontSize',lfs)

%----------------------- fcp -------------------------%
figure(i); i=i+1;
clf 
box  on
hold on
for l=1:Ns
plot(tv,fcp(:,l)*1e2,'--','Color',css(l))
end
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
set(gca,'YLim',ylmfc)
title('Pacific Sediments')
xlabel('Time (y)')
ylabel('CaCO_3 (wt %)')
Hl=legend(dstr);
set(Hl,'FontSize',lfs)

if(ftys)
%----------------------- fct -------------------------%
figure(i); i=i+1;
clf 
box  on
hold on
for l=1:Ns
plot(tv,fct(:,l)*1e2,'--','Color',css(l))
end
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
set(gca,'YLim',ylmfc)
title('Tethys Sediments')
xlabel('Time (y)')
ylabel('CaCO_3 (wt %)')
Hl=legend(dstr);
set(Hl,'FontSize',lfs)

end

%----------------------- ccd -------------------------%
figure(i); i=i+1;
clf 
box  on
hold on
plot(tv,ccda,'- ','Color',ocncol(1))
plot(tv,ccdi,'--','Color',ocncol(2))
plot(tv,ccdp,'-.','Color',ocncol(3))
if(ftys)
plot(tv,ccdt,'-.','Color',ocncol(4))
end
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
set(gca,'YDir','r')
xlabel('Time (y)')
ylabel('CCD (km)')
Hl=legend(ccdstr);
set(Hl,'FontSize',lfs)

if(f13csd)
%----------------------- f13csed -----------------------%
figure(i); i=i+1;
dcca = (fcca./fca/0.011-1.)*1e3;
dcci = (fcci./fci/0.011-1.)*1e3;
dccp = (fccp./fcp/0.011-1.)*1e3;
if(ftys)
dcct = (fcct./fct/0.011-1.)*1e3;
end
clf 
box  on
hold on
for l=1:Ns
plot(tv,dcca(:,l),'- ','Color',ocncol(1))
plot(tv,dcci(:,l),'--','Color',ocncol(2))
plot(tv,dccp(:,l),'-.','Color',ocncol(3))
if(ftys)
plot(tv,dcct(:,l),'-.','Color',ocncol(4))
end
end
hold off
set(gca,'FontSize',fs)
set(gca,'XLim',xlm)
title('C13 Sediments')
xlabel('Time (y)')
ylabel('')
Hl=legend(dstr);
set(Hl,'FontSize',lfs)
end % 13csd

end % fsed

end % plotflag

return
