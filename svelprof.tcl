puts "enter: psf dcd binwidth nequi t(fs) freq fname"
package require pbctools

proc svelprof { psf dcd binwidth nequi t freq {fname "profile.dat"} } {
    
    mol load psf $psf.psf
    mol addfile $dcd.dcd waitfor all
    set nframe [molinfo top get numframes]
    set nprod [expr double($nframe-$nequi)]; #double is very imp 
    set dt [expr 1.0e-6*$t*$freq];  # time in ns	
    set pi 3.1415926535897931;
    
    pbc unwrap -all

    set selcnt [atomselect top "type CA"];
    set mmcnt [measure minmax $selcnt];
    set cmcnt [measure center $selcnt];

    set radius [expr 0.5*([lindex $mmcnt 1 0]-[lindex $mmcnt 0 0])] 
    #set radius [expr [lindex $mmcnt 1 1]-[lindex $mmcnt 0 1]] 
    set zlen   [expr [lindex $mmcnt 1 2]-[lindex $mmcnt 0 2]] 

    set nbin [expr int($radius/$binwidth)+1]
    set totvol [expr 0.001*$pi*$radius*$radius*$zlen]; # nm^3

    for {set k 0} {$k < $nbin} {incr k} {
    set dens($k) 0
    set svel($k) 0
    set volume($k) [expr 0.001*(2*$k+1)*$pi*$binwidth*$binwidth*$zlen]; # nm^3
    }

    for {set i $nequi} {$i < $nframe} {incr i} {
 
          molinfo top set frame $i	
          set nextframe [expr $i+1]
   
          for {set k 0} {$k < $nbin} {incr k} {
          set den($k) 0
          set vel($k) 0
          }



          set water [atomselect top "type OT"]    
          foreach atom [$water get index] {
          set selatom [atomselect top "index $atom"]
          set pos [lindex [$selatom get {x y z} ] 0]
          $selatom frame $nextframe
          set nextpos [lindex [$selatom get {x y z} ] 0]
          $selatom delete

          set x [lindex $pos 0]
          set y [lindex $pos 1]
          set z1 [lindex $pos 2]
          set z2 [lindex $nextpos 2]
          set dist [expr pow (($x*$x + $y*$y), 0.5)];    # assuming the pore axis along z direction 
          set bin  [expr int($dist/$binwidth)]
          incr den($bin)
          set vel($bin) [expr $vel($bin)+0.1*($z2-$z1)];  # distance in nm 
          }
         $water delete

          for {set k 0} {$k < $nbin} {incr k} {
          set dens($k) [expr $dens($k)+$den($k)]
          set svel($k) [expr $svel($k)+$vel($k)]
          }
   }
 
    set tot_1 0
    set tot_2 0
    
    for {set k 0} {$k < $nbin} {incr k} {
    set  tot_1 [expr $tot_1 + $dens($k)]
    set  tot_2 [expr $tot_2 + $volume($k)]
    }
    puts " total water: [expr $tot_1/$nprod]  density: [expr $tot_1/($nprod*$totvol)]  bin_volume: $tot_2 tot_vol: $totvol"

    for {set k 0} {$k < $nbin} {incr k} {
    if { $dens($k)>0 }  {
    set dens($k) [expr $dens($k)/($nprod)]
    set svel($k) [expr $svel($k)/($dt*$nprod*$dens($k))] 
    set dens($k) [expr $dens($k)/($volume($k))] }
    }

    set fp [open $fname w]
    for {set k 0} {$k < $nbin} {incr k} {
    puts $fp "[format %.4f [expr $k*$binwidth+0.5*$binwidth]]  \t [format %.4f $dens($k)] \t [format %.4f $svel($k)]"
    puts "[format %.4f [expr $k*$binwidth+0.5*$binwidth]]  \t [format %.4f $dens($k)]  \t [format %.4f $svel($k)] "
    }
    close $fp

}




