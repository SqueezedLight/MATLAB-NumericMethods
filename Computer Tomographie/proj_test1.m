function profilwertv=proj_test1(xiv,phi)

% Dieses Programm gehoert zu  "ct2.m"  (DeVries "L und Kreis" MIT VERSCHIEDENEN HOEHEN).

% Argument  xiv  kann ein Vektor sein!!

  nxi=length(xiv);
  profilwertv=zeros(1,nxi);

  for i=1:nxi
    xi=xiv(i);

%  Test: Kreis mit Mittelpunkt (xz,yz) und Radius rad:
            xz=0.5;
            yz=0.5;
            rad=0.25;

            if(phi==pi/2.0)
              phi=0.99*phi;
            end
            if(phi==pi)
              phi=0.99*phi;
            end

            cosin=cos(phi);
            sinus=sin(phi);
            alpha=xi*sinus + xi*cosin*cosin/sinus - yz;
            a=1.0+(cosin/sinus)^2;
            b=-2.0*(xz+alpha*cosin/sinus);
            c=xz^2+alpha^2-rad^2;
            root_arg=b^2-4.0*a*c;
            if(root_arg >= 0.0)   
              x1=(-b-sqrt(root_arg))/(2*a);
              x2=(-b+sqrt(root_arg))/(2*a);
              y1=xi*sinus - cosin/sinus*(x1-xi*cosin);
              y2=xi*sinus - cosin/sinus*(x2-xi*cosin);
              profilwert=sqrt((x1-x2)^2+(y1-y2)^2);
            else
              profilwert=0.0;
            end

%  L-Stueck, aus zwei Balken bestehend:
          cosin=cos(phi);
          sinus=sin(phi);
          for ijk=1:2
            if(ijk == 1)   
              xl=-0.6;
              xr=-0.4;
              yu=-0.5;
              yo=0.7;
              height=0.75;
            else 
              xl=-0.4;
              xr=0.3;
              yu=-0.5;
              yo=-0.3;
              height=0.5;
            end

              yl=xi/sinus - cosin/sinus * xl;
              yr=xi/sinus - cosin/sinus * xr;
% Es gibt sieben verschiedene Schnittpunkt-Arten: 
              if(yl>yo & yr>yo) 
                proj=0.0;
              end
              if(yl>yo & yr<=yo & yr>=yu)
                  xo=sinus/cosin*(xi/sinus-yo);
                  proj=sqrt((xo-xr)^2+(yo-yr)^2);
              end
              if(yl<=yo & yl>=yu & yr<=yo & yr>=yu) 
                proj=sqrt((xr-xl)^2+(yr-yl)^2);
              end
              if(yl>yo & yr<yu) 
                xo=sinus/cosin*(xi/sinus-yo);
                xu=sinus/cosin*(xi/sinus-yu);
                proj=sqrt((xo-xu)^2+(yo-yu)^2);
              end
              if(yl<yu & yr>yo)
                xo=sinus/cosin*(xi/sinus-yo);
                xu=sinus/cosin*(xi/sinus-yu);
                proj=sqrt((xo-xu)^2+(yo-yu)^2);
              end
              if(yr<=yu & yl<=yo & yl>=yu) 
                xu=sinus/cosin*(xi/sinus-yu);
                proj=sqrt((xl-xu)^2+(yl-yu)^2);
              end
              if(yl<yu & yr<yu) 
                proj=0.0;
              end
              if(yr>yo & yl<=yo & yl>=yu) 
                xo=sinus/cosin*(xi/sinus-yo);
                proj=sqrt((xl-xo)^2+(yl-yo)^2);
              end
              if(yl<yu & yr<=yo & yr>=yu) 
                xu=sinus/cosin*(xi/sinus-yu);
                proj=sqrt((xu-xr)^2+(yu-yr)^2);
              end
            profilwert=profilwert+height*proj;
          end  

    profilwertv(i)=profilwert;
  end
