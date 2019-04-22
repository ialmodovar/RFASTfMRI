####***********************************************************************
##
## @file: biweight.R
##
## Compute robust estimate of the scale parameter using biweight
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
##
## Authors:
## Ranjan Maitra, 
## Department of Statistics
## Iowa State University
## maitra@iastate.edu
##
## Israel A Almodovar-Rivera, 
## Department of Biostatistics and Epidemiology
## Graduate School of Public Health
## University of Puerto Rico
## Medical Science Campus
## israel.almodovar@upr.edu 
##***********************************************************************

biweight.scale.est <- function(x, center = median(x), w = 6){
    y <- x-center
    n <- length(y)
    stilde <- mad(x)
    u <- y/w/stilde
    num <- sum(y^2*(1-u^2)^4*(u^2 < 1))
    den <- sum((1-u^2)*(1-5*u^2)*(u^2 < 1))
    sqrt(n * num)/abs(den)
}
