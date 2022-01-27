def index_shift(i=None,c="",s="1"):
    """
    Returns a string that describes the gridpoint
    index, allowing for shifts.

    Expected output:

    index_shift():
    "i0_i1_i2"

    index_shift(0,"p"):
    "i0p1_i1_i2"

    index_shift(2,"m","half"):
    "i0_i1_i2mhalf"
    """
    if c == "":
        return "i0_i1_i2"
    else:
        if i == 0:
            return "i0"+c+s+"_i1_i2"
        elif i == 1:
            return "i0_i1"+c+s+"_i2"
        else:
            return "i0_i1_i2"+c+s

def gf_read(gf,i,c=""):
    """
    Returns a string to read a gridfunction from
    main memory.

    Expected output:

    gf_read("<gf>",2,c="p"):
    'const REAL <gf>_i0_i1_i2p1 = <gf>[i0_i1_i2p1];\n'
    """
    index = index_shift(i,c)
    string = "const REAL "+gf+"_"+index+" = "+gf+"["+index+"];\n"
    return string

def read_in_metric(indent="  "):
    """
    Returns a string to read in metric gridfunctions
    from main memory according to the necessities of
    this C function.
    
    Expected output of read_in_metric():
    
    const REAL gammaDD00_i0_i1_i2 = gammaDD00[i0_i1_i2];
    const REAL gammaDD00_i0p1_i1_i2 = gammaDD00[i0p1_i1_i2];
    const REAL gammaDD00_i0m1_i1_i2 = gammaDD00[i0m1_i1_i2];

    const REAL gammaDD11_i0_i1_i2 = gammaDD11[i0_i1_i2];
    const REAL gammaDD11_i0_i1p1_i2 = gammaDD11[i0_i1p1_i2];
    const REAL gammaDD11_i0_i1m1_i2 = gammaDD11[i0_i1m1_i2];

    const REAL gammaDD22_i0_i1_i2 = gammaDD22[i0_i1_i2];
    const REAL gammaDD22_i0_i1_i2p1 = gammaDD22[i0_i1_i2p1];
    const REAL gammaDD22_i0_i1_i2m1 = gammaDD22[i0_i1_i2m1];
    """
    string = ""
    for i in range(3):
        gf = "gammaDD"+str(i)+str(i)
        string += indent+gf_read(gf,i)
        string += indent+gf_read(gf,i,"p")
        string += indent+gf_read(gf,i,"m")+"\n"
    return string

def read_in_opacities(indent="  "):
    """
    Returns a string to read in opacity gridfunctions
    from main memory according to the necessities of
    this C function.

    Expected output of read_in_opacities() is six copies
    of the following lines:

    const REAL <gf>_i0_i1_i2 = <gf>[i0_i1_i2];
    const REAL <gf>_i0p1_i1_i2 = <gf>[i0p1_i1_i2];
    const REAL <gf>_i0m1_i1_i2 = <gf>[i0m1_i1_i2];
    const REAL <gf>_i0_i1p1_i2 = <gf>[i0_i1p1_i2];
    const REAL <gf>_i0_i1m1_i2 = <gf>[i0_i1m1_i2];
    const REAL <gf>_i0_i1_i2p1 = <gf>[i0_i1_i2p1];
    const REAL <gf>_i0_i1_i2m1 = <gf>[i0_i1_i2m1];

    where <gf> is replaced by kappa_0_nue, kappa_1_nue,
    kappa_0_anue, kappa_1_anue, kappa_0_nux, and kappa_1_nux.
    """
    string = ""
    gfs = ["kappa_0_nue","kappa_1_nue","kappa_0_anue","kappa_1_anue","kappa_0_nux","kappa_1_nux"]
    for gf in gfs:
        string += indent+gf_read(gf,0)
        for i in range(3):
            string += indent+gf_read(gf,i,"p")
            string += indent+gf_read(gf,i,"m")
        string += "\n"
    return string

def average_gf(gf,i,c):
    """
    Returns a string that computes the average of
    a given gridfunction. Example outputs:
    
    average_gf("<gf>",0,"p")
    'const REAL <gf>_i0phalf_i1_i2 = 0.5*(<gf>_i0_i1_i2 + <gf>_i0p1_i1_i2);\n'
    
    average_gf("<gf>",2,"m")
    'const REAL <gf>_i0_i1_i2mhalf = 0.5*(<gf>_i0_i1_i2 + <gf>_i0_i1_i2m1);\n'
    """
    return "const REAL "+gf+"_"+index_shift(i,c,"half")+" = 0.5*("+gf+"_"+index_shift(i)+" + "+gf+"_"+index_shift(i,c)+");\n"

def average_metric(indent="  "):
    """
    Returns a string that computes the average metric.

    Expected output of average_metric():

    const REAL gammaDD00_i0phalf_i1_i2 = 0.5*(gammaDD00_i0_i1_i2 + gammaDD00_i0p1_i1_i2);
    const REAL gammaDD00_i0mhalf_i1_i2 = 0.5*(gammaDD00_i0_i1_i2 + gammaDD00_i0m1_i1_i2);

    const REAL gammaDD11_i0_i1phalf_i2 = 0.5*(gammaDD11_i0_i1_i2 + gammaDD11_i0_i1p1_i2);
    const REAL gammaDD11_i0_i1mhalf_i2 = 0.5*(gammaDD11_i0_i1_i2 + gammaDD11_i0_i1m1_i2);

    const REAL gammaDD22_i0_i1_i2phalf = 0.5*(gammaDD22_i0_i1_i2 + gammaDD22_i0_i1_i2p1);
    const REAL gammaDD22_i0_i1_i2mhalf = 0.5*(gammaDD22_i0_i1_i2 + gammaDD22_i0_i1_i2m1);
    """
    string = ""
    for i in range(3):
        gf = "gammaDD"+str(i)+str(i)
        string += indent+average_gf(gf,i,"p")
        string += indent+average_gf(gf,i,"m")+"\n"
    return string

def average_opacities(indent = "  "):
    """
    Returns a string that computes the average opacity.

    Expected output of average_opacities() is n copies of the following
    lines:

    const REAL <gf>_i0phalf_i1_i2 = 0.5*(<gf>_i0_i1_i2 + <gf>_i0p1_i1_i2);
    const REAL <gf>_i0mhalf_i1_i2 = 0.5*(<gf>_i0_i1_i2 + <gf>_i0m1_i1_i2);
    const REAL <gf>_i0_i1phalf_i2 = 0.5*(<gf>_i0_i1_i2 + <gf>_i0_i1p1_i2);
    const REAL <gf>_i0_i1mhalf_i2 = 0.5*(<gf>_i0_i1_i2 + <gf>_i0_i1m1_i2);
    const REAL <gf>_i0_i1_i2phalf = 0.5*(<gf>_i0_i1_i2 + <gf>_i0_i1_i2p1);
    const REAL <gf>_i0_i1_i2mhalf = 0.5*(<gf>_i0_i1_i2 + <gf>_i0_i1_i2m1);

    where <gf> is replaced by each of the n gridfunctions in
    the gfs list defined below.
    """
    string = ""
    gfs = ["kappa_0_nue","kappa_1_nue","kappa_0_anue","kappa_1_anue","kappa_0_nux","kappa_1_nux"]
    for gf in gfs:
        for i in range(3):
            string += indent+average_gf(gf,i,"p")
            string += indent+average_gf(gf,i,"m")
        string += "\n"
    return string

def compute_all_ds(indent="  "):
    """
    Returns a string that sets ds according to

    ds = sqrt( gamma_{ii}dxx^{i}dxx^{i} ).

    Expected output of compute_all_ds():

    const REAL ds_i0phalf_i1_i2 = sqrt(dxx0*dxx0*gammaDD00_i0phalf_i1_i2);
    const REAL ds_i0mhalf_i1_i2 = sqrt(dxx0*dxx0*gammaDD00_i0mhalf_i1_i2);
    const REAL ds_i0_i1phalf_i2 = sqrt(dxx1*dxx1*gammaDD11_i0_i1phalf_i2);
    const REAL ds_i0_i1mhalf_i2 = sqrt(dxx1*dxx1*gammaDD11_i0_i1mhalf_i2);
    const REAL ds_i0_i1_i2phalf = sqrt(dxx2*dxx2*gammaDD22_i0_i1_i2phalf);
    const REAL ds_i0_i1_i2mhalf = sqrt(dxx2*dxx2*gammaDD22_i0_i1_i2mhalf);
    """
    string = ""
    for i in range(3):
        for c in ["p","m"]:
            index = index_shift(i,c,"half")
            string += indent+"const REAL "
            string += "ds_"+index+" = sqrt(dxx"+str(i)+"*dxx"+str(i)+"*gammaDD"+str(i)+str(i)+"_"+index+");\n"
    return string

def compute_tau(indent="  "):
    """
    Returns a string to compute optical depth.

    Expected output of compute_tau() is 6 copies of the following
    lines:

    const REAL tau_0_nue_i0p1_i1_i2 = tau_0_nue[i0p1_i1_i2] + ds_i0phalf_i1_i2*kappa_0_nue_i0phalf_i1_i2;
    const REAL tau_0_nue_i0m1_i1_i2 = tau_0_nue[i0m1_i1_i2] + ds_i0mhalf_i1_i2*kappa_0_nue_i0mhalf_i1_i2;
    const REAL tau_0_nue_i0_i1p1_i2 = tau_0_nue[i0_i1p1_i2] + ds_i0_i1phalf_i2*kappa_0_nue_i0_i1phalf_i2;
    const REAL tau_0_nue_i0_i1m1_i2 = tau_0_nue[i0_i1m1_i2] + ds_i0_i1mhalf_i2*kappa_0_nue_i0_i1mhalf_i2;
    const REAL tau_0_nue_i0_i1_i2p1 = tau_0_nue[i0_i1_i2p1] + ds_i0_i1_i2phalf*kappa_0_nue_i0_i1_i2phalf;
    const REAL tau_0_nue_i0_i1_i2m1 = tau_0_nue[i0_i1_i2m1] + ds_i0_i1_i2mhalf*kappa_0_nue_i0_i1_i2mhalf;

    where we compute the optical depth for number and energy of
    electron neutrinos, antineutrinos, and heavy lepton neutrinos
    and antineutrinos.
    """
    string = ""
    gfs_tau   = ["tau_0_nue"  ,"tau_1_nue"  ,"tau_0_anue"  ,"tau_1_anue"  ,"tau_0_nux"  ,"tau_1_nux"]
    gfs_kappa = ["kappa_0_nue","kappa_1_nue","kappa_0_anue","kappa_1_anue","kappa_0_nux","kappa_1_nux"]
    for k in range(len(gfs_tau)):
        tau   = gfs_tau[k]
        kappa = gfs_kappa[k]
        for i in range(3):
            for c in ["p","m"]:
                index_1 = index_shift(i,c)
                index_2 = index_shift(i,c,"half")
                string += indent+"const REAL "+tau+"_"+index_1+" = "+tau+"["+index_1+"] + ds_"+index_2+"*"+kappa+"_"+index_2+";\n"
        string += "\n"
    return string

def compute_new_tau(indent="  "):
    """
    Compute optical depth using a "path of least resistance"
    algorithm. This means that we look at the neighboring
    values of the optical depth and pick the smallest one.
    """
    gfs = ["tau_0_nue","tau_1_nue","tau_0_anue","tau_1_anue","tau_0_nux","tau_1_nux"]
    string = ""
    for gf in gfs:
        string += indent+"const REAL new_"+gf+"_"+index_shift()+" = "
        string += "MIN(MIN(MIN(MIN(MIN("
        string += gf+"_"+index_shift(0,"p")+","
        string += gf+"_"+index_shift(0,"m")+"),"
        string += gf+"_"+index_shift(1,"p")+"),"
        string += gf+"_"+index_shift(1,"m")+"),"
        string += gf+"_"+index_shift(2,"p")+"),"
        string += gf+"_"+index_shift(2,"m")+");\n"
    return string

def write_tau_to_main_memory(indent="  "):
    """
    Returns a string to write the result of the optical
    depth to main memory.

    Expected output of write_tau_to_main_memory():

    tau_0_nue[i0_i1_i2] = new_tau_0_nue_i0_i1_i2;
    tau_1_nue[i0_i1_i2] = new_tau_1_nue_i0_i1_i2;
    tau_0_anue[i0_i1_i2] = new_tau_0_anue_i0_i1_i2;
    tau_1_anue[i0_i1_i2] = new_tau_1_anue_i0_i1_i2;
    tau_0_nux[i0_i1_i2] = new_tau_0_nux_i0_i1_i2;
    tau_1_nux[i0_i1_i2] = new_tau_1_nux_i0_i1_i2;
    """
    gfs = ["tau_0_nue","tau_1_nue","tau_0_anue","tau_1_anue","tau_0_nux","tau_1_nux"]
    string = ""
    for gf in gfs:
        index = index_shift()
        string += indent+gf+"["+index+"] = new_"+gf+"_"+index+";\n"
    return string

def write_optical_depth_Cfunc_body(indent="  "):
    """
    This function writes the body of the C function
    NRPyLeakage_compute_optical_depths().
    """
    body  = indent+"// Step 1: Set gridpoint indices\n"
    body += indent+"const int i0_i1_i2   = IDX3D(i0  ,i1,i2  );\n"
    body += indent+"const int i0p1_i1_i2 = IDX3D(i0+1,i1,i2  );\n"
    body += indent+"const int i0m1_i1_i2 = IDX3D(i0-1,i1,i2  );\n"
    body += indent+"const int i0_i1p1_i2 = IDX3D(i0,i1+1,i2  );\n"
    body += indent+"const int i0_i1m1_i2 = IDX3D(i0,i1-1,i2  );\n"
    body += indent+"const int i0_i1_i2p1 = IDX3D(i0,i1  ,i2+1);\n"
    body += indent+"const int i0_i1_i2m1 = IDX3D(i0,i1  ,i2-1);\n\n"
    body += indent+"// Step 2: Read in metric gfs from main memory\n"+read_in_metric(indent)+"\n"
    body += indent+"// Step 3: Read in opacity gfs from main memory\n"+read_in_opacities(indent)+"\n"
    body += indent+"// Step 4: Compute metric at cell faces\n"+average_metric(indent)+"\n"
    body += indent+"// Step 5: Compute ds^{i} = sqrt(gamma_{ii}dx^{i}dx^{i})\n"+compute_all_ds(indent)+"\n"
    body += indent+"// Step 6: Compute opacities at cell faces\n"+average_opacities(indent)+"\n"
    body += indent+"// Step 7: Compute optical depth at neighboring points\n"+compute_tau(indent)+"\n"
    body += indent+"// Step 8: Select path of least resistance\n"+compute_new_tau(indent)+"\n"
    body += indent+" // Step 9: Write results to main memory\n"+write_tau_to_main_memory(indent)
    return body
