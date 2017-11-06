/********************************************************************/
/* This is a gauss procedure to compute Moon and Perron (2004),     */
/* "TESTING FOR A UNIT ROOT IN PANELS WITH DYNAMIC FACTORS",        */
/* Journal of Econometrics 122, 81-126.                             */
/*                                                                  */
/* This code is heavily based on the Matlab code by                 */
/* Moon & Perron 2004, last updated June 25 2004.                   */
/*								                                    */
/* {testvec}=PU_MP04(y,m,case);	                                    */
/*                                                                  */
/*  <<INPUT>>                                                       */
/*  y: (T+1 x N) matrix                                             */
/*  m: number of factors (estimated beforehand or known)            */
/*  case:  1 or 2: fixed effects no trend, or 3: deterministic trend*/
/*                                                                  */
/*  <<OUTPUT>>                                                      */
/*  testvec: a (3 x 1) vector, (m,tstara,tstarb)'                   */
/*          tstara and tstarb are defined in p.92 for case 1 & 2.   */
/*                                                                  */
/*  NB:This code uses Andrews and Monahan (1992) bandwidth estimator*/
/* NOTEs: Default setting is to print the results. It can be        */
/*       suppressed by choosing "outpt=0" below                     */
/*  Takashi Yamagata, 6 March 06                                    */
/********************************************************************/


proc(1)=PU_MP04(x_mat,m,case);
local outpt,n,t,al,p,y_mat,y_1_mat,z,x1,qx,ztilde,ztildet,ztildel,rhohat,yhat;
local values,vectors,zqbeta,zqbetal,resid,sigmas,sigmasq,omega;
local lambda1,matj1,omegas,omegasq,matj2,lambda,phi4;
local zq,zql,rhostar,hatrho,nt,tstara,tstarb,testvec;

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
outpt=1;/* 1: reports ourput, 0:supress the output */
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
outwidth 256;

z=x_mat;
n=cols(z);
t=rows(z)-1;

    if case==1 or case==2;
    
    al = 0;    /* a1=0 implies no trend*/
    p = -1;    /* no deterministic component estimated  */ 
    
    elseif case == 3;
        al = 1;  /* linear trend*/
        p = 1;   /* trend estimated*/
    endif;    
    
              /*remove deterministic components*/
        
        if p == 0; 
            x1 = ones(t+1,1);
            qx = eye(t+1)-x1*x1'/(t+1);
        elseif p == 1;
            x1 = ones(t+1,1)~seqa(1,1,t+1);
            qx = eye(t+1)-x1*inv(x1'*x1)*x1';
        elseif p == -1;
            qx = eye(t+1);  
        endif;

        ztilde = qx*z;        
        ztildet = ztilde[2:t+1,.];
        ztildel = ztilde[1:t,.];        
        
        rhohat = tr(ztildel*ztildet')/tr(ztildel*ztildel');
        yhat = ztildet-rhohat*ztildel;       
        
        /*extract factors using normalization with smaller dimension*/
    
        if n <= t;
            {values,vectors} = eighv(yhat'*yhat);              
        else;        
            {values,vectors} = eighv(yhat*yhat');     
        endif;
        
/* number of factors should be given as "m"
*/       
            /*compute FM estimator*/

            zqbeta = ortho(ztildet,yhat,vectors,m);      /* project ztilde onto orthogonal space */
            zqbetal = ortho(ztildel,yhat,vectors,m);     
            
            resid = zqbeta-rhohat*zqbetal;
        
            /*compute variance of errors*/
        
            sigmas = sumc(abs(resid).^2)/t;
            sigmasq = sumc(sigmas)/n;
                        
            /*compute long-run covariances*/
            
            /*(bandwidth selection based on Andrews and Monahan 1992 */  
            omega=zeros(n,1);lambda1=zeros(n,1);                   
            for j (1,n,1);
                {matj1,matj2} = lrcov(resid[.,j],1,999);                
                omega[j]=matj1;
                lambda1[j]=matj2;
            endfor;            
            omegas = omega;
            omegasq = meanc(omegas);
            lambda = 0.5*(omegasq - sigmasq);
            phi4 = (omegas'*omegas)/n;            
                
            zq = zqbetal*zqbeta';        
            zql = zqbetal*zqbetal';        
            
            /* autoregressive parameter - rhostar with estimated variances */        
            
            if p == -1 or p==0; 
                rhostar = (tr(zq)-n*t*meanc(diag(lambda)))/tr(zql); 
            elseif p==1;
                rhostar = (tr(zq)+0.5*n*t*meanc(sigmas))/tr(zql);
            endif;

            /*compute test statistics*/
            
            NT = sqrt(n)*t;
            if p == -1 or p==0;
                tstara = NT*(rhostar-1)/sqrt(2*phi4/omegasq^2);
                tstarb = NT*(rhostar-1)*sqrt(tr(zql)/(n*t^2))*sqrt(omegasq/phi4);    
            elseif p == 1;
                tstara = NT*(rhostar-1)/sqrt((15/4)*phi4/omegasq^2);
                tstarb = NT*(rhostar-1)*sqrt(tr(zql)/(n*t^2))*sqrt(4*omegasq/phi4);
            endif;  
testvec=(m|tstara|tstarb);

/**/
if outpt==1;
"";
"MOON AND PERRON (2004) PANEL UNIT ROOT TEST RESULTS";
if case==1 or case==2;
"In the case without trend";
elseif case==3;
"In the case with trend";
endif;
format /rd 4,0;
"";
"N=";;n;
"Number of data points used:";;t+1;;"(T = ";;t;;")";
"Given number of factors:";;m;
"";
format /rd 10,3;
"tstar_a";;tstara;
"tstar_b";;tstarb;
"------------------------------------------------";
endif;
retp(testvec);
endp;     




/***********************************************************************/
/*************   PROCEDURES ********************************************/
/***********************************************************************/

proc(1)=tr(x);
local tr;
tr=sumc(diag(x));
retp(tr);
endp;

proc(2) = LRcov(x,wght,mtd);
local t,n,xt,xl,sigma,a,v,d,i,xpw,omega,b,resid,s_resid,s4,wght1;
local s_deux,s_zero,delta,l,lambda,j,gama,w,qs,dinv;

t=rows(x);n=cols(x);
if n==1;
wght1=wght;
else;
wght1=wght|ones(n-1,1);
endif;
/*prewhitening*/

sigma = x'*x/t;

xt = x[2:t,.];
xl = x[1:t-1,.];

a = xt'*xl*inv(xl'*xl);
{d,v} = eighv(a);
for i (1,n,1);
    if (d[i,i]>0.97) or (d[i,i]<-0.97);
        d[i,i] = 0.97*(2*(d[i,i].>1)-1);
    endif;
endfor;
a = v*d*inv(v);

xpw = xt - xl*a';

omega = xpw'*xpw/t;

/*------------------*/
if mtd==999;
    /* Andrews*/
    xt = xpw[2:t-1,.];
    xl = xpw[1:t-2,.];
    b = diag((xt'*xl)./(xl'*xl));
    Resid = xt-xl*diag(b);        
    s_resid = Resid'*Resid/t;
    s4 = diag(s_resid)^2;
    
    s_deux = (4*((b)^2)*s4/((1-b)^8));
    s_zero = ( s4 ./ ((1-b).^4)); 
    delta = 1.3221*((wght1'*(s_deux))/(wght1'*(s_zero))).^(0.2);
    
endif;
/*---------------*/
L = delta*t^(0.2);

lambda = 0;
for j (1,t-2,1);
    gama = xpw[j+1:t-1,.]'*xpw[1:t-j-1,.]/t;
    v=j/L;
    w=6*pi*v/5;        
    qs=((25./(12*((pi*v).^2))).*(((sin(w))./w)-cos(w)));
    omega = omega + qs*(gama + gama');
    lambda = lambda +qs*gama;
    
endfor;

/* recolor */

dinv = inv(eye(n)-a);
omega = dinv*omega*dinv';
lambda = dinv*lambda*dinv' + sigma*a'*dinv';

retp (omega,lambda);
endp;


proc(1) = ortho(y,yhat,v,k);
local t,n,loadings,factors,betahat,norm,qbeta,x;

t=rows(y);n=cols(y);

if n <= t;
    loadings=sqrt(n)*v[.,n-k+1:n];
    factors=yhat*loadings/n;
    betahat=(yhat'*yhat)*loadings/(n*t);
else;
    factors=sqrt(t)*v[.,t-k+1:t];
    loadings=yhat'*factors/t;
    norm=loadings'*loadings/n;
    betahat=loadings*chol(norm);
endif;

qbeta=eye(n)-betahat*inv(betahat'*betahat)*betahat';

x = y*qbeta;

retp(x);
endp;
