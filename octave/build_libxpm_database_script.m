#!/usr/bin/octave

arg_list = argv();

if(length(arg_list) != 7)
	disp("check_last_inputs requires 7 parameters");
	return;
end

fs_new = str2num(arg_list{1});
channel_power_new = str2num(arg_list{2});
D_new = str2num(arg_list{3});
alphaDB_new = str2num(arg_list{4});
gam_new = str2num(arg_list{5});
res_Disp_new = str2num(arg_list{6});
HalfWindow_new = str2num(arg_list{7});

build_libxpm_database(fs_new,channel_power_new,D_new,alphaDB_new,gam_new,res_Disp_new,HalfWindow_new);
