using FITSIO.Libcfitsio;
using OIFITS;
using Dates;
#=
function get_extname(filein)
    num=fits_get_num_hdus(filein)
    extnames=Array{String}(num)
    for i=1:num
      hdu_type=fits_movabs_hdu(filein,i)
      if hdu_type ==:binary_table
        extnames[i]=fits_read_keyword(filein,"EXTNAME")[1]
      else
        extnames[i]="NONE"
      end
    end
    return extnames
  end

function find_hdu_from_extname(extname_array,extname_name)
   indx=find(extname_array->extname_array==extname_name,extname_array)
  return indx
end


function make_header_cfitsio(header,fileout)
  for i=1:length(header)
    fits_write_key(fileout,keys(header)[i],values(header)[i],get_comment(keys(header)[i]));
  end
end

function copy_oifits(fileout,filein)
    fts = FITS(filein)
    copy_oi_header(fileout,exten,fts)
    length(fts)       #number of hdu
    for i=2:length
        hdr = FITSIO.read_header(fts[i]])
end

function copycols(fileout,hduin)
    colnames=FITSIO.colnames(hduin)
    for i=1:length(colnames)
        col=FITSIO.read(hduin,colnames(i))
        fits_write_col(fileout, i, 1, 1, col;
    end
end
=#

function copy_oi_header(fileout,hduin) #filenone in this case is in FITSIO not cfitsio
  #= Make header. Based on John Young's code =#
  #Move to primary HDU
  if fits_get_num_hdus(fileout) == 0
    fits_create_img(fileout,Int16,zeros(Int64,0)); #Create primary HDU if it doesn't exist
  else
    fits_movabs_hdu(fileout,1);  #Move to primary HDU
  end
  header=FITSIO.read_header(hduin);
  nkeys=length(header)
    for i = 4:nkeys #first four are created by create_img
        entry = FITSIO.read_key(hduin,i);
        if typeof(entry[2]) != Nothing
            fits_write_key(fileout,entry[1],entry[2],entry[3]);
        else
        fits_write_key(fileout,entry[1],"",entry[3]);
        end
    end
 end

function copy_oi_array(fileout,hduin,arrays)

  #=
  Writing the OI_WAVELENGTH binary table
  based on cfitsio and John Young's oifitslib
  =#
  header=FITSIO.read_header(hduin);
  tfields = 5;
  ttype= ["TEL_NAME", "STA_NAME", "STA_INDEX", "DIAMETER", "STAXYZ","FOV", "FOVTYPE"];
  tform = ["16A", "16A", "I", "E", "3D", "D", "16A"];
  tunit = ["\0", "\0", "\0", "m", "m", "arcsec", ""];
  extname = header["EXTNAME"];
  extver = header["EXTVER"]
  revision = header["OI_REVN"]; #1 or 2 depending on which format XYZ
  arrname=header["ARRNAME"]
  frame=header["FRAME"]
  arrayx=header["ARRAYX"]
  arrayy=header["ARRAYY"]
  arrayz=header["ARRAYZ"]
  if revision == 2
    tfields = 7;
      coldefs = [(ttype[1],tform[1],tunit[1]),(ttype[2],tform[2],tunit[2]),(ttype[3],tform[3],tunit[3]),
            (ttype[4],tform[4],tunit[4]),(ttype[5],tform[5],tunit[5]),(ttype[6],tform[6],tunit[6]),(ttype[7],tform[7],tunit[7])];
          #= OFITS2 needs to add ,(ttype[6],tform[6],tunit[6]),(ttype[7],tform[7],tunit[7])];=#
          fits_create_binary_tbl(fileout,0,coldefs,extname);
          fits_write_key(fileout,"EXTNAME",extname,"name of the binary extension")
          fits_write_key(fileout,"EXTEND",extver,"Extension version");
          fits_write_key(fileout,"OI_REVN",revision,"Revision number of the table definition");
          fits_write_key(fileout,"ARRNAME",arrname,"Array Name");
          fits_write_key(fileout,"FRAME",frame,"Coordinate Frame");
          fits_write_key(fileout,"ARRAYX",arrayx,"[m] Array center x coordinate");
          fits_write_key(fileout,"ARRAYY",arrayy,"[m] Array center y coordinate");
          fits_write_key(fileout,"ARRAYZ",arrayz,"[m] Array center z coordinate");
          fits_write_col(fileout, 1, 1, 1, convert(Array{String},arrays[1]));
          fits_write_col(fileout, 2, 1, 1, convert(Array{String},arrays[2]));
          fits_write_col(fileout, 3, 1, 1, convert(Array{Int16},arrays[3]));
          fits_write_col(fileout, 4, 1, 1, convert(Array{Float32},arrays[4]));
          fits_write_col(fileout, 5, 1, 1, convert(Array{Float64,2},arrays[5]));
          fits_write_col(fileout, 6, 1, 1, convert(Array{Float64,1},arrays[6]));
          fits_write_col(fileout, 7, 1, 1, convert(Array{String,1},arrays[7]));
          #= OIFITS2 STUFF
          int revision = 2, irow;
          fits_write_key(fptr, TSTRING, "ARRNAOI_VIS2 ME", array.arrname,
            "Array name", pStatus);
  fits_write_key(fptr, TSTRING, "FRAME", array.frame,
                 "Coordinate frame", pStatus);
  fits_write_key(fptr, TDOUBLE, "ARRAYX", &array.arrayx,
                 "Array centre x coordinate", pStatus);
  fits_write_key_unit(fptr, "ARRAYX", "m", pStatus);
  fits_write_key(fptr, TDOUBLE, "ARRAYY", &array.arrayy,
                 "Array centre y coordinate", pStatus);
  fits_write_key_unit(fptr, "ARRAYY", "m", pStatus);
  fits_write_key(fptr, TDOUBLE, "ARRAYZ", &array.arrayz,
                 "Array centre z coordinate", pStatus);
  fits_write_key_unit(fptr, "ARRAYZ", "m", pStatus);
  fits_write_key(fptr, TINT, "EXTVER", &extver,
"ID number of this OI_ARRAY", pStatus);=#
#  irows=size(OIFITS.select(oifilein,"OI_ARRAY"))[1]

  #OIFITS2 needs additional columns which might be more difficult as OIFITS doesn't have the get procedures for those yet
    #if error status ==1 return 0; this is an error check I might put in later
    else
        coldefs = [(ttype[1],tform[1],tunit[1]),(ttype[2],tform[2],tunit[2]),(ttype[3],tform[3],tunit[3]),
              (ttype[4],tform[4],tunit[4]),(ttype[5],tform[5],tunit[5])];
            #= OFITS2 needs to add ,(ttype[6],tform[6],tunit[6]),(ttype[7],tform[7],tunit[7])];=#
            fits_create_binary_tbl(fileout,0,coldefs,extname);
            fits_write_key(fileout,"EXTNAME",extname,"name of the binary extension")
            fits_write_key(fileout,"EXTEND",extver,"Extension version");
            fits_write_key(fileout,"OI_REVN",revision,"Revision number of the table definition");
            fits_write_key(fileout,"ARRNAME",arrname,"Array Name");
            fits_write_key(fileout,"FRAME",frame,"Coordinate Frame");
            fits_write_key(fileout,"ARRAYX",arrayx,"[m] Array center x coordinate");
            fits_write_key(fileout,"ARRAYY",arrayy,"[m] Array center y coordinate");
            fits_write_key(fileout,"ARRAYZ",arrayz,"[m] Array center z coordinate");
            fits_write_col(fileout, 1, 1, 1, convert(Array{String},arrays[1]));
            fits_write_col(fileout, 2, 1, 1, convert(Array{String},arrays[2]));
            fits_write_col(fileout, 3, 1, 1, convert(Array{Int16},arrays[3]));
            fits_write_col(fileout, 4, 1, 1, convert(Array{Float32},arrays[4]));
            fits_write_col(fileout, 5, 1, 1, convert(Array{Float64,2},arrays[5]));
        end
 end


function copy_oi_target(fileout,hduin,arrays)
  #=
  /**
   * Write OI_TARGET fits binary table
   */=#
  header=FITSIO.read_header(hduin);
  tfields = 17;
  ttype = ["TARGET_ID", "TARGET", "RAEP0", "DECEP0", "EQUINOX",
                   "RA_ERR", "DEC_ERR", "SYSVEL", "VELTYP",
                   "VELDEF", "PMRA", "PMDEC", "PMRA_ERR", "PMDEC_ERR",
                   "PARALLAX", "PARA_ERR", "SPECTYP"];
  tform = ["I", "16A", "D", "D", "E",
                   "D", "D", "D", "8A",
                   "8A", "D", "D", "D", "D",
                   "E", "E", "16A"];
  tunit = ["\0", "\0", "deg", "deg", "yr",
                   "deg", "deg", "m/s", "\0",
                   "\0", "deg/yr", "deg/yr", "deg/yr", "deg/yr",
                   "deg", "deg", "\0"];
  extname = header["EXTNAME"];
  revision =header["OI_REVN"] ;
  coldefs =[(ttype[1],tform[1],tunit[1]),(ttype[2],tform[2],tunit[2]),(ttype[3],tform[3],tunit[3]),
          (ttype[4],tform[4],tunit[4]),(ttype[5],tform[5],tunit[5]),(ttype[6],tform[6],tunit[6]),
          (ttype[7],tform[7],tunit[7]),(ttype[8],tform[8],tunit[8]),(ttype[9],tform[9],tunit[9]),
          (ttype[10],tform[10],tunit[10]),(ttype[11],tform[11],tunit[11]),(ttype[12],tform[12],tunit[12]),
          (ttype[13],tform[13],tunit[13]),(ttype[14],tform[14],tunit[14]),(ttype[15],tform[15],tunit[15]),
          (ttype[16],tform[16],tunit[16]),(ttype[17],tform[17],tunit[17])];
  fits_create_binary_tbl(fileout,0,coldefs,extname);
  fits_write_key(fileout,"EXTNAME",extname,"name of the binary extension")
  fits_write_key(fileout,"OI_REVN",revision,"Revision number of the table definition");
  fits_write_col(fileout, 1, 1, 1,  convert(Array{Int16},arrays[1]));
  fits_write_col(fileout, 2, 1, 1,  convert(Array{String},arrays[2]));
  fits_write_col(fileout, 3, 1, 1,  convert(Array{Float64},arrays[3]));
  fits_write_col(fileout, 4, 1, 1,  convert(Array{Float64},arrays[4]));
  fits_write_col(fileout, 5, 1, 1,  convert(Array{Float32},arrays[5]));
  fits_write_col(fileout, 6, 1, 1,  convert(Array{Float64},arrays[6]));
  fits_write_col(fileout, 7, 1, 1,  convert(Array{Float64},arrays[7]));
  fits_write_col(fileout, 8, 1, 1,  convert(Array{Float64},arrays[8]));
  fits_write_col(fileout, 9, 1, 1,  convert(Array{String},arrays[9]));
  fits_write_col(fileout, 10, 1, 1, convert(Array{String},arrays[10]));
  fits_write_col(fileout, 11, 1, 1, convert(Array{Float64},arrays[11]));
  fits_write_col(fileout, 12, 1, 1, convert(Array{Float64},arrays[12]));
  fits_write_col(fileout, 13, 1, 1, convert(Array{Float64},arrays[13]));
  fits_write_col(fileout, 14, 1, 1, convert(Array{Float64},arrays[14]));
  fits_write_col(fileout, 15, 1, 1, convert(Array{Float32},arrays[15]));
  fits_write_col(fileout, 16, 1, 1, convert(Array{Float32},arrays[16]));
  fits_write_col(fileout, 17, 1, 1, convert(Array{String},arrays[17]));
  #fits_write_key(fileout,"EXTEND",extver,"Extension version");
  #= OIFITS2 Stuff
  int revision = 2, irow
  OIFITS2 needs additional columns which might be more difficult as OIFITS doesn't have the get procedures for those yet =#
  #if error status ==1 return 0; this is an error check I might put in later
 end

function copy_oi_wavelength(fileout,hduin,arrays)
  #=
  Writing the OI_WAVELENGTH binary table
  based on cfitsio and John Young's oifitslib
  =#
  header=FITSIO.read_header(hduin);
  tfields = 2;
  ttype= ["EFF_WAVE", "EFF_BAND"];
  tform= ["1E", "1E"];
  tunit = ["m", "m"];
  extname = header["EXTNAME"];
  extver = header["EXTVER"]
  revision = header["OI_REVN"]; #1 or 2 depending on which format XYZ
  insname=header["INSNAME"]

  coldefs =[(ttype[1],tform[1],tunit[1]),(ttype[2],tform[2],tunit[2])];
  fits_create_binary_tbl(fileout,0,coldefs,extname);
  fits_write_key(fileout,"EXTNAME",extname,"name of the binary extension")
  fits_write_key(fileout,"EXTEND",extver,"Extension version")
 fits_write_key(fileout,"INSNAME",insname,"Name of detector, for cross-referencing");
  fits_write_key(fileout,"OI_REVN",revision,"Revision number of the table definition");

  #=OIFITS2 STUFF
  CHECK REVISION? if (wave.revision != revision) { printf("WARNING! wave.revision != %d on entry to %s. ""Writing revision %d table\n", revision, function, revision);}
  FITSIO.Libcfitsio.fits_write_key(fout,"OI_REVN", revision,"Revision number of the table definition");
  FITSIO.Libcfitsio.fits_write_key(fout, "INSNAME", wave.insname,"Detector name");
  FITSIO.Libcfitsio.fits_write_key(fout,"EXTVER",extver,"ID number of this OI_WAVELENGTH"); #Normally 1 it seems
    #if error status ==1 return 0; this is an error check I might put in later
  #nwave=size(OIFITS.select(oifilein,"OI_WAVELENGTH"))[1]#
=#

  fits_write_col(fileout, 1, 1, 1, convert(Array{Float32},arrays[1]));
  fits_write_col(fileout, 2, 1, 1, convert(Array{Float32},arrays[2]));
 end


function copy_oi_vis2(fileout,hduin,arrays)
#=
/**
 * Write OI_VIS2 fits binary table
 *
 * @param fptr     see cfitsio documentation
 * @param vis2     data struct, see exchange.h
 * @param extver   value for EXTVER keyword
 * @param pStatus  pointer to status variable
 *
 * @return On error, returns non-zero cfitsio error code, and sets *pStatus
 */
=#
  header=FITSIO.read_header(hduin);
  tfields = 10;  # mandatory columns */
  zerotime = 0.0;
  nw=size(arrays[5])[1]
  ttype = ["TARGET_ID", "TIME", "MJD", "INT_TIME",
           "VIS2DATA", "VIS2ERR", "UCOORD", "VCOORD",
           "STA_INDEX", "FLAG"];
  #tform = ["1I", "1D", "1D", "1D",
    #       "8D", "8D", "1D", "1D",
    #        "2I", "8L"];
    tform = ["1I", "1D", "1D", "1D",
             "$(nw)D", "$(nw)D", "1D", "1D",
             "2I", "$(nw)L"];
  tunit = ["\0", "s", "day", "s",
            "\0", "\0", "m", "m",
            "\0", "\0"];
  date_obs  = string(Dates.today())
  date_obs = header["DATE-OBS"]
  arrname = header["ARRNAME"]
  extname = header["EXTNAME"];
  extver = header["EXTVER"]
  revision = header["OI_REVN"]; #1 or 2 depending on which format XYZ
  insname=header["INSNAME"]
  coldefs =[(ttype[1],tform[1],tunit[1]),(ttype[2],tform[2],tunit[2]),(ttype[3],tform[3],tunit[3]),
          (ttype[4],tform[4],tunit[4]),(ttype[5],tform[5],tunit[5]),(ttype[6],tform[6],tunit[6]),
          (ttype[7],tform[7],tunit[7]),(ttype[8],tform[8],tunit[8]),(ttype[9],tform[9],tunit[9]),
          (ttype[10],tform[10],tunit[10])];
  fits_create_binary_tbl(fileout,0,coldefs,extname);
 fits_write_key(fileout,"EXTNAME",extname,"name of the binary extension")
  fits_write_key(fileout,"EXTEND",extver,"Extension version");
  fits_write_key(fileout,"OI_REVN",revision,"Revision number of the table definition");
  fits_write_key(fileout,"DATE-OBS",date_obs,"UTC start date of observations");
  fits_write_key(fileout,"ARRNAME",arrname,"(optional) Identifies corresponding OI_Array");
  fits_write_key(fileout,"INSNAME",insname,"Identifies corresponding OI_WAVELENGTH");
  fits_write_key(fileout,"EXTEND",extver,"Extension version");
  fits_write_col(fileout, 1, 1, 1,  convert(Array{Int16},arrays[1]));
  fits_write_col(fileout, 2, 1, 1,  convert(Array{Float64},arrays[2]));
  fits_write_col(fileout, 3, 1, 1,  convert(Array{Float64},arrays[3]));
  fits_write_col(fileout, 4, 1, 1,  convert(Array{Float64},arrays[4]));
  fits_write_col(fileout, 5, 1, 1,  convert(Array{Float64},arrays[5]));
  fits_write_col(fileout, 6, 1, 1,  convert(Array{Float64},arrays[6]));
  fits_write_col(fileout, 7, 1, 1,  convert(Array{Float64},arrays[7]));
  fits_write_col(fileout, 8, 1, 1,  convert(Array{Float64},arrays[8]));
  fits_write_col(fileout, 9, 1, 1,  convert(Array{Int16},arrays[9]));
  fits_write_col(fileout, 10, 1, 1, convert(Array{Bool,2},arrays[10]));
 end

function copy_oi_t3(fileout,hduin,arrays)
#=
/**
 * Write OI_T3 fits binary tableDATE
 *
 * @param fptr     see cfitsio documentation
 * @param vis2     data struct, see exchange.h
 * @param extver   value for EXTVER keyword
 * @param pStatus  pointer to status variable
 *
 * @return On error, returns non-zero cfitsio error code, and sets *pStatus
 */
=#
  header=FITSIO.read_header(hduin);
  nw=size(arrays[5])[1]
  tfields = 14;  # mandatory columns */
  zerotime = 0.0;
  ttype = ["TARGET_ID", "TIME", "MJD", "INT_TIME",
           "T3AMP", "T3AMPERR", "T3PHI","T3PHIERR",
            "U1COORD", "V1COORD","U2COORD", "V2COORD",
           "STA_INDEX", "FLAG"];
  #tform = ["1I", "1D", "1D", "1D",
    #       "8D", "8D", "8D", "8D",
    #       "1D","1D","1D","1D","3I", "8L"];

    tform = ["1I", "1D", "1D", "1D",
            "$(nw)D", "$(nw)D", "$(nw)D", "$(nw)D",
            "1D","1D","1D","1D","3I", "$(nw)L"];
  tunit = ["\0", "s", "day", "s",
           "\0", "\0", "deg", "deg",
           "m", "m", "m", "m",
           "\0", "\0"];

  date_obs = header["DATE-OBS"]
  arrname = header["ARRNAME"]
  extname = header["EXTNAME"];
  extver = header["EXTVER"]
  revision = header["OI_REVN"]; #1 or 2 depending on which format XYZ
  insname=header["INSNAME"]

  coldefs =[(ttype[1],tform[1],tunit[1]),(ttype[2],tform[2],tunit[2]),(ttype[3],tform[3],tunit[3]),
          (ttype[4],tform[4],tunit[4]),(ttype[5],tform[5],tunit[5]),(ttype[6],tform[6],tunit[6]),
          (ttype[7],tform[7],tunit[7]),(ttype[8],tform[8],tunit[8]),(ttype[9],tform[9],tunit[9]),
          (ttype[10],tform[10],tunit[10]),(ttype[11],tform[11],tunit[11]),(ttype[12],tform[12],tunit[12]),
          (ttype[13],tform[13],tunit[13]),(ttype[14],tform[14],tunit[14])];
  fits_create_binary_tbl(fileout,0,coldefs,extname);
    fits_write_key(fileout,"EXTNAME",extname,"name of the binary extension")
  fits_write_key(fileout,"OI_REVN",revision,"Revision number of the table definition");
  fits_write_key(fileout,"DATE-OBS",date_obs,"UTC start date of observations");
  fits_write_key(fileout,"ARRNAME",arrname,"(optional) Identifies corresponding OI_Array");
  fits_write_key(fileout,"INSNAME",insname,"Identifies corresponding OI_WAVELENGTH");
  fits_write_key(fileout,"EXTEND",extver,"Extension version");

  #/* Write mandatory columns */
  fits_write_col(fileout, 1, 1, 1,  convert(Array{Int16},arrays[1]));
  fits_write_col(fileout, 2, 1, 1,  convert(Array{Float64},arrays[2]));
  fits_write_col(fileout, 3, 1, 1,  convert(Array{Float64},arrays[3]));
  fits_write_col(fileout, 4, 1, 1,  convert(Array{Float64},arrays[4]));
  fits_write_col(fileout, 5, 1, 1,  convert(Array{Float64},arrays[5]));
  fits_write_col(fileout, 6, 1, 1,  convert(Array{Float64},arrays[6]));
  fits_write_col(fileout, 7, 1, 1,  convert(Array{Float64},arrays[7]));
  fits_write_col(fileout, 8, 1, 1,  convert(Array{Float64},arrays[8]));
  fits_write_col(fileout, 9, 1, 1,  convert(Array{Float64},arrays[9]));
  fits_write_col(fileout, 10, 1, 1,  convert(Array{Float64},arrays[10]));
  fits_write_col(fileout, 11, 1, 1,  convert(Array{Float64},arrays[11]));
  fits_write_col(fileout, 12, 1, 1,  convert(Array{Float64},arrays[12]))
  fits_write_col(fileout, 13, 1, 1,  convert(Array{Int16},arrays[13]));
  fits_write_col(fileout, 14, 1, 1, convert(Array{Bool,2},arrays[14]));
 end
