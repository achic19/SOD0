import os
import osmnx as ox

project_crs = 'epsg:3857'
import warnings
import geopandas as gpd
import ast

warnings.filterwarnings(action='ignore')
pjr_loc = os.getcwd()


## Choose locations and create folders if necessary
class Preprocessing:
    def __init__(self, place):
        self.place = place
        self.round_about = None  # for all the elements of roundabout type
        self.is_junction = True  # True if our dataset include the column  'junction' after downloading from OSM
        self.data_folder = None  # all the new shpfiles will be stored in here

    def create_folder(self, is_test=False):
        """
        This function run over the place and create the necessary folders
        :param place:
        :param is_test:
        :return:
        """

        if is_test:
            self.data_folder = f'{os.path.dirname(os.path.dirname(pjr_loc))}' \
                               f'/places/{self.place.replace(",", "_").replace(" ", "_")}_test'
            os.makedirs(f'{self.data_folder}/delete_2_nodes', exist_ok=True)
            os.makedirs(f'{self.data_folder}/split_tp_intersection', exist_ok=True)
            os.makedirs(f'{self.data_folder}/simplification', exist_ok=True)
        else:
            self.data_folder = f'{os.path.dirname(pjr_loc)}/places/{self.place.replace(",", "_").replace(" ", "_")}'
        print(f'data folder:{self.data_folder}')
        os.makedirs(f'{self.data_folder}', exist_ok=True)
        return self.data_folder

    def first_filtering(self):
        """
        First polylines filtering based on some criteria along
        :return:
        """
        ## Filter out polylines and calculate angles
        my_gdf = gpd.read_file(f'{self.data_folder}/osm_data.gpkg',
                               layer='edges')  # Identify roundabout elements, if any exist, and store them in a
        # separate DataFrame.
        if self.place == 'Tel Aviv':
            my_gdf.rename(columns={'name:en': 'name'}, inplace=True)

        self.is_junction = True if 'junction' in my_gdf.columns else False
        if self.is_junction:
            self.round_about = my_gdf[my_gdf['junction'].isin(['roundabout', 'circular'])]
            my_gdf = my_gdf[~((my_gdf['junction'] == 'roundabout') | (my_gdf['junction'] == 'circular'))]
        to_remove = my_gdf[~((my_gdf['highway'] == 'motorway') | (my_gdf['highway'] == 'trunk') | (
                my_gdf['highway'] == 'motorway_link') | (my_gdf['highway'] == 'motorway_link') | (
                                     my_gdf['highway'] == 'trunk_link'))]
        # Eliminate polylines that lack a name and calculate angles ranging from 0 to 180 degrees based on the
        # bearing field.
        df_pro = to_remove.to_crs(project_crs).dropna(subset=['name'])
        df_pro = df_pro[df_pro['name'] != '']
        df_pro['angle'] = df_pro['bearing'].apply(lambda x: x if x < 180 else x - 180)
        df_pro['length'] = df_pro.length

        # Function to convert valid list strings to lists
        def convert_to_list(s):
            try:
                return ast.literal_eval(s)[0]
            except (ValueError, SyntaxError, TypeError):
                return s  # Return the original string if conversion fails

        # Apply the function to the DataFrame column so polylines with several street names will return the first name.
        df_pro['name'] = df_pro['name'].apply(convert_to_list)
        df_pro['highway'] = df_pro['highway'].apply(convert_to_list)

        # Remove streets that are not connected to any other streets
        # vars
        str_name = 'name_left'
        con_str_name = 'name_right'
        df_analysis = df_pro.copy()
        s_join_analysis = gpd.sjoin(df_analysis, df_pro)
        s_join_analysis2 = s_join_analysis[s_join_analysis[str_name] != s_join_analysis[con_str_name]]
        not_connected = set(df_pro['name']) - set(s_join_analysis2.reset_index()[str_name])
        df_pro = df_pro[~df_pro['name'].isin(not_connected)]
        return df_pro



