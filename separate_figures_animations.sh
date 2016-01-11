#!/bin/bash

DATADIR=~/Documents/DSP_UROP/all_data

# create directories
mkdir -p $DATADIR/animations/command_healthy
mkdir -p $DATADIR/animations/command_impaired
mkdir -p $DATADIR/animations/copy_healthy
mkdir -p $DATADIR/animations/copy_impaired

mkdir -p $DATADIR/figures/command_healthy
mkdir -p $DATADIR/figures/command_impaired
mkdir -p $DATADIR/figures/copy_healthy
mkdir -p $DATADIR/figures/copy_impaired

mkdir -p $DATADIR/figures/command_healthy_velocity_acceleration
mkdir -p $DATADIR/figures/command_impaired_velocity_acceleration
mkdir -p $DATADIR/figures/copy_healthy_velocity_acceleration
mkdir -p $DATADIR/figures/copy_impaired_velocity_acceleration

mkdir -p $DATADIR/figures/command_healthy_spectra
mkdir -p $DATADIR/figures/command_impaired_spectra
mkdir -p $DATADIR/figures/copy_healthy_spectra
mkdir -p $DATADIR/figures/copy_impaired_spectra

mkdir -p $DATADIR/figures/command_healthy_pressure
mkdir -p $DATADIR/figures/command_impaired_pressure
mkdir -p $DATADIR/figures/copy_healthy_pressure
mkdir -p $DATADIR/figures/copy_impaired_pressure

# transfer files
#cp $DATADIR/figs_raw/YDU*/command*mp4 $DATADIR/animations/command_healthy
#cp $DATADIR/figs_raw/CIN*/command*mp4 $DATADIR/animations/command_impaired
#cp $DATADIR/figs_raw/YDU*/copy*mp4 $DATADIR/animations/copy_healthy
#cp $DATADIR/figs_raw/CIN*/copy*mp4 $DATADIR/animations/copy_impaired

#cp $DATADIR/figs_raw/YDU*/xy_command*png $DATADIR/figures/command_healthy
#cp $DATADIR/figs_raw/CIN*/xy_command*png $DATADIR/figures/command_impaired
#cp $DATADIR/figs_raw/YDU*/xy_copy*png $DATADIR/figures/copy_healthy
#cp $DATADIR/figs_raw/CIN*/xy_copy*png $DATADIR/figures/copy_impaired

#cp $DATADIR/figs_raw/YDU*/vxt_axt_command*png $DATADIR/figures/command_healthy_velocity_acceleration
#cp $DATADIR/figs_raw/CIN*/vxt_axt_command*png $DATADIR/figures/command_impaired_velocity_acceleration
#cp $DATADIR/figs_raw/YDU*/vxt_axt_copy*png $DATADIR/figures/copy_healthy_velocity_acceleration
#cp $DATADIR/figs_raw/CIN*/vxt_axt_copy*png $DATADIR/figures/copy_impaired_velocity_acceleration

#cp $DATADIR/figs_raw/YDU*/dft*_command*png $DATADIR/figures/command_healthy_spectra
#cp $DATADIR/figs_raw/CIN*/dft*_command*png $DATADIR/figures/command_impaired_spectra
#cp $DATADIR/figs_raw/YDU*/dft*_copy*png $DATADIR/figures/copy_healthy_spectra
#cp $DATADIR/figs_raw/CIN*/dft*_copy*png $DATADIR/figures/copy_impaired_spectra

cp $DATADIR/figs_raw/YDU*/pt_command*png $DATADIR/figures/command_healthy_pressure
cp $DATADIR/figs_raw/CIN*/pt_command*png $DATADIR/figures/command_impaired_pressure
cp $DATADIR/figs_raw/YDU*/pt_copy*png $DATADIR/figures/copy_healthy_pressure
cp $DATADIR/figs_raw/CIN*/pt_copy*png $DATADIR/figures/copy_impaired_pressure

echo Done
