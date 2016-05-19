#!/bin/bash
 
# svn-up
#
# purpose:  
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail.com
# web:      http://ocefpaf.tiddlyspot.com/
# date:     22-Dec-09
# modified: 22-Dec-09
#
# obs: download all toolboxes
#

# check dependencies
commds=( svn )
for commd in ${commds[@]}; do
  if [[ -z $( type -p ${commd} ) ]]; then
    echo -e "${commd} -- NOT INSTALLED !"
    exit 1
  fi
done

tools=( "http://mirone.googlecode.com/svn/trunk/ mirone"
        "https://mexcdf.svn.sourceforge.net/svnroot/mexcdf/mexnc/trunk mexnc"
        "https://mexcdf.svn.sourceforge.net/svnroot/mexcdf/snctools/branches/nc4/ snctools"
        "https://mexcdf.svn.sourceforge.net/svnroot/mexcdf/netcdf_toolbox/trunk netcdf_toolbox"
        "https://www.myroms.org/svn/om/matlab/roms_wilkin matlab_wilkin"
        "https://www.myroms.org/svn/src/matlab matlab_arango"
        "https://svn1.hosted-projects.com/cmgsoft/m_cmg/trunk/seagrid  seagrid"
        "https://svn1.hosted-projects.com/cmgsoft/m_cmg/trunk/cmgtool/CMGmfiles/ CMGmfiles"
        "https://svn1.hosted-projects.com/cmgsoft/m_cmg/trunk/cmglib/ cmglib"
        "http://svn1.hosted-projects.com/cmgsoft/m_contrib/trunk/rslice/ rslice"
        "http://svn1.hosted-projects.com/cmgsoft/m_cmg/trunk/njTools/release/beta/2.0 njToolbox-2.0"
        "http://svn1.hosted-projects.com/cmgsoft/m_cmg/trunk/RPSstuff RPSstuff"
        "http://svn1.hosted-projects.com/cmgsoft/m_contrib/trunk/seawater seawater"
        "http://svn1.hosted-projects.com/cmgsoft/m_contrib/trunk/tidal_ellipse tidal_ellipse"
        "http://svn1.hosted-projects.com/cmgsoft/m_contrib/trunk/t_tide t_tide"
        "http://svn1.hosted-projects.com/cmgsoft/m_contrib/trunk/air_sea air_sea"
        "http://svn1.hosted-projects.com/cmgsoft/m_contrib/trunk/rdadcp rdadcp"
        "http://imos-toolbox.googlecode.com/svn/trunk/ imos-toolbox"
        "http://wafo.googlecode.com/svn/trunk/ wafo"
      )

for tool in "${tools[@]}"; do
    svn co $tool
done

# http://svn1.hosted-projects.com/cmgsoft/m_contrib/trunk/
# ABS_Toolbox AquatecAquascat1000 LP_Bathymetry MuellerADCP_PD02MAT Roms_tools SequoiaLISST100 SequoiaLISSTST darren_stuff diwasp ezyfit gridfitdir herbers lsge mdd mexcdf_2008b_plus misc specstuff wafo whgrp_adcp

# https://svn1.hosted-projects.com/cmgsoft/m_cmg/trunk/
# ABSS_tbx AQDlib AWAC LISST MMstuff abslib adcirc_tides adcp_tbx advlib cf cmglib  cmgtool dolly fvcom hydratools meta_tools midas ndbc njTools omviz roms_crs roms_jcw roms_rps seagauge_tools sonarlib stpete swan_crs timeplt tri wave_tbx whfc_to_ge