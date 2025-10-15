# -*- coding: utf-8 -*-

from os.path import basename, join
from time import sleep

import arcpy
from arcpy.sa import ExtractByMask, ExtractMultiValuesToPoints

FIELD_CALC_DEFINITIONS = """
def absolute_diff_wsels(wsel1, wsel2):
    if wsel1 > 0 and wsel2 > 0: #ignores nodata values
        return round(abs(wsel1 - wsel2), 1)
    else:
        return -9999
def valid_diff_wsels(abs_diff, tolerance):
    if abs_diff == -9999:
        return "Missed"
    elif abs_diff <= tolerance:
        return f"Passed {tolerance}"
    else:
        return f"Failed {tolerance}"
def get_wsel_value(wsel1, wsel2):
    return wsel1 + wsel2
"""
DATATYPE_SHP = [
    "DEShapefile",
    "DEDbaseTable",
    "DEFeatureClass",
    "GPFeatureLayer",
    "GPLayer",
    "DETable", 
    "GPTableView",
    "DEDbaseTable",
    "DEArcInfoTable",
    "GPString",
]
DATATYPE_TIF = [
    "DEMosaicDataset",
    "DERasterBand",
    "DERasterCatalog",
    "GPRasterCatalogLayer",
    "GPRasterDataLayer",
    "DERasterDataset",
    "GPRasterLayer",
]
DATATYPE_DIR = ["DEFolder","DEWorkspace","GPString"]
DATATYPE_FLD = "Field"

FBS2DTest_HELP_DESCRIPTION = ( 
    "Working Standard ID 128: For floodplains mapped from 2D models, "
    + "evaluation lines and BFE lines on the FIRM must match modeled water "
    + "surface elevations and must be plotted at intervals sufficient to interpolate "
    + "accurate BFEs in between BFE or evaluation lines. If this is not possible, "
    + "separate Flood Profiles for significant flow paths and/or FIS Report inserts "
    + "must also be created.\n"
    + "According to FEMA’s Mapping Base Flood Elevations on FIRMs Guidance "
    + "(November 2023): Where BFEs and evaluation lines are the primary source "
    + "of water surface elevation information on the FIRM, the mapping partner "
    + "should confirm that the lines placed are sufficient to linearly interpolate "
    + "the BFE at any point in the floodplain within 0.1 foot of the true value from "
    + "the model water surface elevation grid.\n"
    + "STEP 1 | After placing BFEs and evaluation lines in accordance with FEMA’s "
    + "Guidelines & Standards, contractors should compare a water surface elevation "
    + "(WSEL) grid generated from the mapped evaluation lines and BFEs to the raw WSEL "
    + "grid output from the model to determine if the BFE and evaluation line spacing "
    + "is adequate (i.e., that interpolated water surface elevations fall within 0.1 "
    + "foot of the modeled results).\n"
    + "STEP 2 | The mapped WSEL grid can be compared to the modeled WSEL grid by "
    + "extracting WSEL values to 'Test Points' and calculating the difference between "
    + "the two surfaces. When creating Test Points, it is important to ensure that "
    + "the number or density of test points is sufficient to provide a reliable "
    + "assessment of the entire mapped grid. At a minimum, Test Points must be "
    + "reasonably and uniformly distributed across the mapped floodplain at a density "
    + "of at least 1 point per 2500 square feet (50x50 feet grid), with a visible increase "
    + "in density through critical areas such as urban or developed areas, floodways, "
    + "areas with significant elevation changes, and regions with high flood risk. This "
    + "is achieved through exporting the computation points from the HEC-RAS model.\n"
    + "STEP 3 | Once the difference in the WSEL grids at the test point locations "
    + "has been calculated, the points are categorized as passing or failing "
    + "depending on if the difference is within the 0.1-foot tolerance. The passing "
    + "percentage should be calculated on a stream-by-stream basis. "
    + "Detailed studies should pass with a reliability of 95% whereas Approximate "
    + "studies should pass with a reliability of 90%."
    + " If the stream passes and the "
    + "community does not have a preference on developing a Flood Profile, no "
    + "further action is needed besides the documentation mentioned in the "
    + "following section.\n"
    + "STEP 4 |If the stream is not passing, the BFEs and evaluation lines should "
    + "be supplemented until the 0.1-foot tolerance is achieved. The points can "
    + "be symbolized by the unique pass/fail values in order to clearly see areas "
    + "within the floodplain that are failing and need additional BFEs. Please "
    + "refer to the Mapping Base Flood Elevations on FIRMs Guidance for further "
    + "details on placing BFEs.\n"
    + "In some situations, it will not be possible to provide the required detail "
    + "using BFE and evaluation lines due to spatial constraints and/or steep "
    + "slope in the water surface profile. In these situations, especially if "
    + "the water surface profile is not linear between mapped BFE or evaluation "
    + "lines, the contractor and CTP should contact the Region 8 Engineer and "
    + "PTS to determine whether Flood Profiles and/or FIS inserts  will be "
    + "needed to provide the required detail for the product user. If a Flood "
    + "Profile or FIS insert is considered inappropriate for the study reach, "
    + "an exception to working standard ID 128 may be considered.\n"
)


# Helpful function to create functional tin and raster from GIS lines
def xs_to_tin_to_raster_by_wtrnm(xs, f_profile, dem, output_dir):
    """
    takes cross sections with WSE values and converts them to a raster
    """
    arcpy.env.snapRaster = dem
    arcpy.env.outputCoordinateSystem = dem
    spatial_ref = arcpy.Describe(dem).spatialReference
    cellsizex = arcpy.GetRasterProperties_management(
        in_raster=dem,
        property_type="CELLSIZEX")

    outtin_name = "SFHA_WSEL"
    outtin = join(output_dir, outtin_name)
    outtif = join(output_dir, "{}.tif".format(outtin_name))

    arcpy.env.workspace = output_dir
    arcpy.ddd.CreateTin(
        out_tin=outtin_name,
        spatial_reference=spatial_ref,
        in_features=[[xs, f_profile, "Hard_Line", "<None>"]],
        constrained_delaunay="DELAUNAY",
    )
    arcpy.ddd.TinRaster(
        in_tin=outtin_name,
        out_raster=outtif,
        data_type="FLOAT",
        method="LINEAR",
        sample_distance="CELLSIZE", # {}".format(cellsizex),
        z_factor=1,
        sample_value=cellsizex)

    return outtif, outtin


def merge_eval_lines(xslines, wsel_field, bfelines, elev_field, out_directory):
    eval_lines = join(out_directory, "S_WSEL_Lines.shp")

    arcpy.management.Merge(
        inputs=[xslines, bfelines],
        output=eval_lines,
        )
    arcpy.management.CalculateField(
        in_table=eval_lines,
        field="GIS_ELEV",
        expression=f"get_wsel_value(!{wsel_field}!,!{elev_field}!)",
        expression_type="PYTHON3",
        code_block=FIELD_CALC_DEFINITIONS,
        field_type="DOUBLE",
    )

    return eval_lines


def sid_128_2d_test( 
        xslines,
        wsel_field,
        bfelines,
        elev_field,
        sfha_01pct,
        ras_wsel_01pct,
        tolerance,
        grid_size,
        out_directory,
        exclusion_polys,
        exclusion_comment_field,
        sfha_stillwater_field,
    ):
    """
    Using FEMA Region 8 2D Approach to test model grids against evaluation lines
    """
    arcpy.AddMessage("Determine environment settings from RAS WSEL grid")
    arcpy.env.outputCoordinateSystem = ras_wsel_01pct
    arcpy.env.snapRaster = ras_wsel_01pct
    arcpy.env.cellSize = ras_wsel_01pct

    remove_datasets = []

    arcpy.AddMessage("Merge the eval XS lines with BFE Lines to map flood fingers")
    if bfelines is not None:
        eval_lines = merge_eval_lines(xslines, wsel_field, bfelines, elev_field, out_directory)
    else:
        eval_lines = join(out_directory, "S_WSEL_Lines.shp")
        arcpy.management.CopyFeatures(in_features=xslines, out_feature_class=eval_lines)
        arcpy.management.CalculateField(
            in_table=eval_lines,
            field="GIS_ELEV",
            expression=f"!{wsel_field}!",
            expression_type="PYTHON3",
            code_block=FIELD_CALC_DEFINITIONS,
            field_type="DOUBLE",
        )

    arcpy.AddMessage("Create Tin and WSEL from GIS to compare with model tif")
    outtif, outtin = xs_to_tin_to_raster_by_wtrnm(
        eval_lines,
        "GIS_ELEV",
        ras_wsel_01pct,
        out_directory
    )

    
    if sfha_stillwater_field:
        arcpy.AddMessage("Stillwater has static field defined, create raster for those areas")
        select_stillwater_areas = join(out_directory, "select_stillwater_areas.shp")
        wc_static = arcpy.AddFieldDelimiters(sfha_01pct, sfha_stillwater_field) + " > -9999"
        arcpy.analysis.Select(
            in_features=sfha_01pct,
            out_feature_class=select_stillwater_areas,
            where_clause = wc_static
        )
        feature_count = int((arcpy.management.GetCount(select_stillwater_areas)).getOutput(0))
        if feature_count > 0:
            select_stillwater_tif = join(out_directory, "select_stilwater.tif")
            arcpy.conversion.PolygonToRaster(
                in_features=select_stillwater_areas, 
                value_field=sfha_stillwater_field, 
                out_rasterdataset=select_stillwater_tif, 
                cell_assignment="CELL_CENTER", 
            )
            outtif_mosaic_name = basename(outtif).replace(".tif", "_mosaic.tif")
            outtif_mosaic_path = join(out_directory, outtif_mosaic_name)
            arcpy.management.MosaicToNewRaster(
                input_rasters=[select_stillwater_tif, outtif], 
                output_location=out_directory, 
                raster_dataset_name_with_extension=outtif_mosaic_name, 
                pixel_type="32_BIT_FLOAT", 
                number_of_bands=1, 
                mosaic_method="FIRST",
            )
            outtif_extract_name = basename(outtif).replace(".tif", "_.tif")
            outtif_extract_path = join(out_directory, outtif_extract_name)
            extract_wsel = ExtractByMask(outtif_mosaic_path, sfha_01pct)
            extract_wsel.save(outtif_extract_path)

            use_outtif = outtif_extract_path
            remove_datasets.extend(
                [select_stillwater_areas, select_stillwater_tif, outtif, outtif_mosaic_path]
            )
        else:
            arcpy.AddMessage("No stillwater areas determined, use tin tif only")
            remove_datasets.append(select_stillwater_areas)
            use_outtif = outtif
            pass
    else: # no stillwater areas found.
        arcpy.AddMessage("No stillwater areas determined, use tin tif only")
        use_outtif = outtif

    arcpy.AddMessage("Grid index Zone A Floodplains to get testing grid.")
    test_mesh_polys = join(out_directory, "S_Grid_Index.shp")
    arcpy.cartography.GridIndexFeatures(
        out_feature_class=test_mesh_polys,
        in_features=sfha_01pct,
        intersect_feature="INTERSECTFEATURE",
        polygon_width=grid_size, #"50 Feet" Default
        polygon_height=grid_size, #"50 Feet" Default
    )

    arcpy.AddMessage("Intersect grid with SFHA to get good test point candidates")
    test_mesh_fp_intersect = join(out_directory, "S_Test_Index_Features.shp")
    arcpy.analysis.Intersect(
        in_features=[test_mesh_polys, sfha_01pct],
        out_feature_class=test_mesh_fp_intersect,
        join_attributes="ONLY_FID",
    )

    arcpy.AddMessage("Get centroids of the intersected grid and mesh.")
    # Increases likelyhood of points being located on model grid and in floodplain
    test_mesh_fp_centroids = join(out_directory, "S_Test_Centroid_Features.shp")
    arcpy.management.FeatureToPoint(
        in_features=test_mesh_fp_intersect,
        out_feature_class=test_mesh_fp_centroids,
        point_location="INSIDE"
    )

    arcpy.AddMessage("Extract values from model grid and gis grid")
    ExtractMultiValuesToPoints(
        in_point_features=test_mesh_fp_centroids, 
        in_rasters=[[ras_wsel_01pct, "RAS_WSEL"],[use_outtif, "SFHA_WSEL"]]
    )
    arcpy.AddMessage("Round extracted values to nearest tenth for simplicity")
    arcpy.management.CalculateField(
        in_table=test_mesh_fp_centroids,
        field="RAS_WSEL",
        expression="round(!RAS_WSEL!, 1)",
        expression_type="PYTHON3",
    )
    arcpy.management.CalculateField(
        in_table=test_mesh_fp_centroids,
        field="SFHA_WSEL",
        expression="round(!SFHA_WSEL!, 1)",
        expression_type="PYTHON3",
    )

    arcpy.AddMessage("Test for tolerance")
    arcpy.management.CalculateField(
        in_table=test_mesh_fp_centroids,
        field="abs_diff",
        expression="absolute_diff_wsels(!RAS_WSEL!, !SFHA_WSEL!)",
        expression_type="PYTHON3",
        code_block=FIELD_CALC_DEFINITIONS,
        field_type="DOUBLE",
    )
    arcpy.management.CalculateField(
        in_table=test_mesh_fp_centroids,
        field="valid_diff",
        expression=f"valid_diff_wsels(!abs_diff!, {tolerance})",
        expression_type="PYTHON3",
        code_block=FIELD_CALC_DEFINITIONS,
        field_type="TEXT",
    )

    out_summary_table = join(out_directory, "T_2D_Summary_Results.dbf")

    if exclusion_polys:
        arcpy.AddMessage("Removing Exclustion Areas")
        test_mesh_fp_centroids_select = join(out_directory, "S_Test_Centroid_Features_Select.shp")
        exclusion_polys_copy = join(out_directory, "S_Exclusion_Areas.shp")
        arcpy.management.CopyFeatures(
            in_features=exclusion_polys,
            out_feature_class=exclusion_polys_copy,
        )
        arcpy.analysis.Erase(
            in_features=test_mesh_fp_centroids,
            erase_features=exclusion_polys_copy,
            out_feature_class=test_mesh_fp_centroids_select,
        )
        arcpy.AddMessage("Create summary statistics for results")
        arcpy.analysis.Statistics(
            in_table=test_mesh_fp_centroids_select,
            out_table=out_summary_table,
            statistics_fields=[["valid_diff", "COUNT"]],
            case_field=["valid_diff"],
        )

    else:
        arcpy.AddMessage("Create summary statistics for results")
        arcpy.analysis.Statistics(
            in_table=test_mesh_fp_centroids,
            out_table=out_summary_table,
            statistics_fields=[["valid_diff", "COUNT"]],
            case_field=["valid_diff"],
        )
        
    arcpy.AddMessage("Analyze results for user") 
    result_messages(
        out_directory,
        out_summary_table,
        tolerance,
        grid_size,
        exclusion_polys,
        exclusion_comment_field,
    )

    remove_datasets.extend([test_mesh_polys, outtin, out_summary_table])

    return remove_datasets


def result_messages(
        out_directory, 
        out_summary_table, 
        tolerance, 
        grid_size, 
        exclusion_polys,
        exclusion_comment_field,
    ):
    """"""
    for row in arcpy.da.SearchCursor(
        in_table=out_summary_table,
        field_names=["valid_diff", "FREQUENCY"]
    ):
        if row[0] == f"Failed {tolerance}":
            failing = int(row[1])
        if row[0] == f"Passed {tolerance}":
            passing = int(row[1])
        if row[0] == "Missed":
            missing = int(row[1])

    total_points = failing + passing + missing
    total_valid_points = failing + passing
    miss_pct = round(float(missing)/total_points, 3) * 100
    valid_pass_pct = round(float(passing)/total_valid_points, 3) * 100
    valid_fail_pct = round(float(failing)/total_valid_points, 3) * 100

    results_txt_file = join(out_directory, 'Results.txt')

    msg = (
        "\n#### FEMA R8 SID 128 Check ####\n \n "
        + "SID 128 2D Test Summary Results\n \n"
        + f"Mesh test size: {grid_size}\n"
        + f"Vertical Test Tolerance: {tolerance}\n \n"
        + f"- Total Points Created: {total_points}\n"
        + f"\t-Missed: {missing} ({miss_pct:0.1f}%)\n \n"
        + f"- Total Valid Points (Not Missed): {total_valid_points}\n"
        + f"\t-Pass: {passing} ({valid_pass_pct:0.1f}%)\n"
        + f"\t-Fail: {failing} ({valid_fail_pct:0.1f}%)\n \n"
    )
    arcpy.AddMessage(msg)
    with open(results_txt_file, "w") as file:
        file.write(msg)
    
    if miss_pct > 10.0:
        msg = "Missed Percentage too High (more than 10%), Failed"
        arcpy.AddMessage(msg)
        with open(results_txt_file, "a") as file:
            file.write(msg)
        return None
    
    if valid_pass_pct >= 95.0:
        msg = "Passed More than 95%\n Detailed with Floodway, Passed"
        arcpy.AddMessage(msg)
        with open(results_txt_file, "a") as file:
            file.write(msg)

    elif 95.0 > valid_pass_pct >= 90.0:
        msg = (
            r"Passed Between 90% and 95%"
            + "\n Detailed with Floodway: Failed\n Detailed without Floodway, Passed"
        )
        arcpy.AddMessage(msg)
        with open(results_txt_file, "a") as file:
            file.write(msg)

    else:
        msg = "Detailed with or without Floodway (less than 90%), Failed"
        arcpy.AddMessage(msg)
        with open(results_txt_file, "a") as file:
            file.write(msg)

    if exclusion_polys:
        exclude_reasons = {}
        for row in arcpy.da.SearchCursor(
            in_table=exclusion_polys,
            field_names=[exclusion_comment_field]
        ):
            if row[0] in exclude_reasons.keys():
                exclude_reasons[row[0]] += 1
            else:
                exclude_reasons[row[0]] = 1
        
        with open(results_txt_file, "a") as file:
            file.write("\n\nExclusion Area Reasons Include:\n")
            for key in exclude_reasons.keys():
                file.write(f"\t{key}: {exclude_reasons[key]} instances\n")

    arcpy.AddMessage(f"Results Summary: {results_txt_file}")

    return None


def input_messages(
        xslines,
        wsel_field,
        bfelines,
        elev_field,
        sfha_01pct,
        ras_wsel_01pct,
        tolerance,
        grid_size,
        out_directory,
        exclusion_polys,
        exclusion_comment_field,
        sfha_stillwater_field,
    ):
    """"""
    arcpy.AddMessage("Input Datasets:\n ")
    arcpy.AddMessage(f"S_XS path = {xslines}")
    arcpy.AddMessage(f"S_XS Field = {wsel_field}")
    arcpy.AddMessage(f"S_BFE path = {bfelines}")
    arcpy.AddMessage(f"S_BFE Field = {elev_field}")
    arcpy.AddMessage(f"S_Fld_Haz_Ar (1% only) path = {sfha_01pct}")
    arcpy.AddMessage(f"S_Fld_Haz_Ar Static_BFE Field = {sfha_stillwater_field}")
    arcpy.AddMessage(f"Model 1% WSEL Grid path = {ras_wsel_01pct}")
    arcpy.AddMessage(f"Tolerance between Model and FEMA datasets = {tolerance}")
    arcpy.AddMessage(f"Testing Grid width/height = {grid_size}")
    arcpy.AddMessage(f"Output directory = {out_directory}")
    arcpy.AddMessage(f"Exclusion Polygons Areas path = {exclusion_polys}")
    arcpy.AddMessage(f"Exclusion Polygons Areas Comment Field = {exclusion_comment_field}\n")
    return None


class Toolbox:
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "FEMA R8 SID 128 Check" 
        self.alias = "fema_r8_sid_128_check"

        # List of tool classes associated with this toolbox
        self.tools = [FBS2DTest]


class FBS2DTest:
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "SID 128 2D Test"
        self.description = FBS2DTest_HELP_DESCRIPTION

    def getParameterInfo(self):
        """Define the tool parameters."""
        xslines = arcpy.Parameter(
                    displayName = "(1) S_XS or S_BFE",
                    name = "xslines",
                    datatype = DATATYPE_SHP,
                    parameterType = "Required",
                    direction = "Input",
            )
        wsel_field = arcpy.Parameter(
                    displayName = "(1) Water Surface Elevation Field",
                    name = "wsel_field",
                    datatype = DATATYPE_FLD,
                    parameterType = "Required",
                    direction = "Input",
            )
        bfelines = arcpy.Parameter(
                    displayName = "(2) Additional S_XS or S_BFE (if applicable)",
                    name = "bfelines",
                    datatype = DATATYPE_SHP,
                    parameterType = "Optional",
                    direction = "Input",
            )
        elev_field = arcpy.Parameter(
                    displayName = "(2) Additional Water Surface Elevation Field (if applicable)",
                    name = "elev_field",
                    datatype = DATATYPE_FLD,
                    parameterType = "Optional",
                    direction = "Input",
            )
        sfha_01pct = arcpy.Parameter(
                    displayName = "(3) S_Fld_Haz_Ar (Zone AE or AH only)",
                    name = "sfha_01pct",
                    datatype = DATATYPE_SHP,
                    parameterType = "Required",
                    direction = "Input",
            )
        sfha_stillwater_field = arcpy.Parameter(
                    displayName = "(3) S_Fld_Haz_Ar STATIC_BFE Field (if applicable)",
                    name = "sfha_stillwater_field",
                    datatype = DATATYPE_FLD,
                    parameterType = "Optional",
                    direction = "Input",
            )
        ras_wsel_01pct = arcpy.Parameter(
                    displayName = r"(4) Model 1% Water Surface Elevation Grid",
                    name = "ras_wsel_01pct",
                    datatype = DATATYPE_TIF,
                    parameterType = "Required",
                    direction = "Input",
            )
        tolerance = arcpy.Parameter(
                    displayName = "(5) Testing Tolerance (units assumed in US Feet)",
                    name = "tolerance",
                    datatype = "GPString",
                    parameterType = "Required",
                    direction = "Input",
            )
        grid_size = arcpy.Parameter(
                    displayName = "(5) Testing Grid (width/height)",
                    name = "grid_size",
                    datatype = "GPLinearUnit",
                    parameterType = "Required",
                    direction = "Input",
            )
        out_directory = arcpy.Parameter(
                    displayName = "(6) Output Directory for Results",
                    name = "out_directory",
                    datatype = DATATYPE_DIR,
                    parameterType = "Required",
                    direction = "Input",
            )
        exclusion_polys = arcpy.Parameter(
                    displayName = "(7) Exclusion Polygon Areas (if applicable)",
                    name = "exclusion_polys",
                    datatype = DATATYPE_SHP,
                    parameterType = "Optional",
                    direction = "Input",
            )
        exclusion_comment_field = arcpy.Parameter(
                    displayName = "(7) Exclusion Polygon Areas Comment Field (if applicable)",
                    name = "exclusion_comment_field",
                    datatype = DATATYPE_FLD,
                    parameterType = "Optional",
                    direction = "Input",
            )
        wsel_field.parameterDependencies = [xslines.name]
        elev_field.parameterDependencies = [bfelines.name]
        sfha_stillwater_field.parameterDependencies = [sfha_01pct.name]
        exclusion_comment_field.parameterDependencies  = [exclusion_polys.name]
        tolerance.value = 0.1 # Vertical Feet Tolerance
        grid_size.value = "50 FEETUS" # U.S. Survey Feet
        params_2d_test = [
            xslines,
            wsel_field,
            bfelines,
            elev_field,
            sfha_01pct,
            sfha_stillwater_field,
            ras_wsel_01pct,
            tolerance,
            grid_size,
            out_directory,
            exclusion_polys,
            exclusion_comment_field,
        ]
        return params_2d_test

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""
        try:
            if arcpy.CheckExtension("3D") != "Available":
                raise Exception
        except Exception:
            return False  # The tool cannot be run

        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False  # The tool cannot be run

        return True

    def updateParameters(self, params_2d_test):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        # Selects Appropriate Fields for S_XS and S_BFE assuming DFIRM schema
        # Otherwise allows user to select field with appropriate data.
        if params_2d_test[0].altered and not params_2d_test[0].hasBeenValidated:
            s_xs = params_2d_test[0].value
            if s_xs:
                field_names = [f.name for f in arcpy.ListFields(s_xs)]
                if "WSEL_REG" in field_names:
                    params_2d_test[1].value = "WSEL_REG"
                elif "ELEV" in field_names:
                    params_2d_test[1].value = "ELEV"
                else:
                    params_2d_test[1].value = ""
        if params_2d_test[2].altered and not params_2d_test[2].hasBeenValidated:
            s_bfe = params_2d_test[2].value
            if s_bfe:
                field_names = [f.name for f in arcpy.ListFields(s_bfe)]
                if "WSEL_REG" in field_names:
                    params_2d_test[3].value = "WSEL_REG"
                elif "ELEV" in field_names:
                    params_2d_test[3].value = "ELEV"
                else:
                    params_2d_test[3].value = ""
        if params_2d_test[4].altered and not params_2d_test[4].hasBeenValidated:
            s_bfe = params_2d_test[4].value
            if s_bfe:
                field_names = [f.name for f in arcpy.ListFields(s_bfe)]
                if "STATIC_BFE" in field_names:
                    params_2d_test[5].value = "STATIC_BFE"
                else:
                    params_2d_test[5].value = ""
        # Selects Comment or comment field for exclusion areas.
        # Otherwise allows user to select field with appropriate data.
        if params_2d_test[10].altered and not params_2d_test[10].hasBeenValidated:
            exclude_zones = params_2d_test[10].value
            if exclude_zones:
                field_names = [f.name for f in arcpy.ListFields(exclude_zones)]
                if "Comment" in field_names:
                    params_2d_test[11].value = "Comment"
                elif "comment" in field_names:
                    params_2d_test[11].value = "comment"
                else:
                    params_2d_test[11].value = ""

        return None

    def updateMessages(self, params_2d_test):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""
        return

    def execute(self, params_2d_test, messages):
        """The source code of the tool."""
        xslines = params_2d_test[0].valueAsText
        wsel_field = params_2d_test[1].valueAsText
        bfelines = params_2d_test[2].valueAsText
        elev_field = params_2d_test[3].valueAsText
        sfha_01pct = params_2d_test[4].valueAsText
        sfha_stillwater_field = params_2d_test[5].valueAsText
        ras_wsel_01pct = params_2d_test[6].valueAsText
        tolerance = float(params_2d_test[7].valueAsText)
        grid_size = params_2d_test[8].valueAsText
        out_directory = params_2d_test[9].valueAsText
        exclusion_polys = params_2d_test[10].valueAsText
        exclusion_comment_field = params_2d_test[11].valueAsText

        input_messages(
            xslines,
            wsel_field,
            bfelines,
            elev_field,
            sfha_01pct,
            ras_wsel_01pct,
            tolerance,
            grid_size,
            out_directory,
            exclusion_polys,
            exclusion_comment_field,
            sfha_stillwater_field,
        )

        arcpy.AddMessage("\n#### Running SID 128 2D Check #####\n\n")
        features = sid_128_2d_test(
            xslines,
            wsel_field,
            bfelines,
            elev_field,
            sfha_01pct,
            ras_wsel_01pct,
            tolerance,
            grid_size,
            out_directory,
            exclusion_polys,
            exclusion_comment_field,
            sfha_stillwater_field,
        )

        for item in features:
            try:
                sleep(1)
                arcpy.management.Delete(item)
            except Exception:
                msg = (
                    f"Unable to remove '{item}'\n"
                    + "Please remove this manually from working directory.\n"
                )
                arcpy.AddWarning(msg)

        return None

    def postExecute(self, params_2d_test):
        """This method takes place after outputs are processed and
        added to the display."""
        return
