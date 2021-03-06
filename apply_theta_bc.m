% ## Copyright (C) 2013 homu
% ## 
% ## This program is free software; you can redistribute it and/or modify
% ## it under the terms of the GNU General Public License as published by
% ## the Free Software Foundation; either version 3 of the License, or
% ## (at your option) any later version.
% ## 
% ## This program is distributed in the hope that it will be useful,
% ## but WITHOUT ANY WARRANTY; without even the implied warranty of
% ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% ## GNU General Public License for more details.
% ## 
% ## You should have received a copy of the GNU General Public License
% ## along with Octave; see the file COPYING.  If not, see
% ## <http://www.gnu.org/licenses/>.

% ## apply_pres_bc

% ## Author: logan-shi <loganshi@sjtu.edu.cn>
% ## Created: 2020-06-07

function [ phi ] = apply_theta_bc (phi)

simple_globals;

% homogeneous BC
%x_low hot
phi(1,:) = T_H*ones(1,ny+1);
%x_high cold
phi(nx+1,:) = T_L*ones(1,ny+1);
%y_low adiabatic
phi(2:nx,1) = phi(2:nx,2);
%y_high adiabatic
phi(2:nx,ny+1) = phi(2:nx,ny);

return
end
