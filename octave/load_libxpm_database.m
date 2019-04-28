#!/usr/bin/octave

function matrix = load_libxpm_database()
  try
    load xpm_database.mat;
  catch
    disp("ERROR: xpm_database.mat file not found or invalid");
    return;
  end
  matrix = xpm_matrix;
endfunction