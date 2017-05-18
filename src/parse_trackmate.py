# -*- coding: utf-8 -*-
"""

This file contains a collection of functions for reading and parsing xml files 
output by the TrackMate ImageJ plugin

Written by Vladislav Belyy (UCSF)

"""
import numpy as np
from xml.dom import minidom

def read_trackmate_file(file):
    """Read a file that was output by the TrackMate plugin.
    
    Args:
        file (str): Full path to TrackMate xml file
        
    Returns:
        parsed_file (Document): TrackMate file parsed with minidom
        
    """
    parsed_file = minidom.parse(file)
    return parsed_file

def get_raw_tracks(parsed_trackmate_file):
    """Get raw tracks from parsed TrackMate file.
    
    Args:
        parsed_trackmate_file (Document): TrackMate file parsed with minidom
    
    Returns:
        tracks_raw (NodeList): DOM node list of all tracks
    """
    tracks_raw = parsed_trackmate_file.getElementsByTagName("Track")
    return tracks_raw

def list_tracks(parsed_trackmate_file):
    """Convert the raw track data to a list of dictionaries.
    
    Args:
        parsed_trackmate_file (Document): TrackMate file parsed with minidom
    
    Returns:
        tracks (list): List of dictionaries containing the attributes of each 
            track.
    """
    tracks = list()
    tracks_raw = get_raw_tracks(parsed_trackmate_file)
    for track_raw in tracks_raw:
        att = dict(track_raw.attributes.items())
        tracks.append(att)
    return tracks
    
def list_spots(parsed_trackmate_file):
    """Convert the raw spot data to a list of dictionaries.
    
    Args:
        parsed_trackmate_file (Document): TrackMate file parsed with minidom
    
    Returns:
        spots (list): List of dictionaries containing the attributes of each 
            spot.
    """
    spots = list()
    spots_by_frame = parsed_trackmate_file.getElementsByTagName("SpotsInFrame")
    for spots_curr_frame in spots_by_frame:
        curr_spots = spots_curr_frame.getElementsByTagName("Spot")
        for spot in curr_spots:
            att = dict(spot.attributes.items())
            spots.append(att)
    return spots       

def pull_property_by_track(prop_name, all_spot_ids, spots):
    """Return any given numerical property of arbitrarily grouped spots.
    
    Args:
        prop_name (str): Any valid key for individal elements in spots.
        all_spot_ids (list): list of numpy arrays (one per group of spots, e.g.
            a single track), each line of the array containing the spot's 
            unique ID and its index in spots
        spots (list): List of dictionaries containing the attributes of each 
            spot.
            
    Returns:
        all_track_property (list): List of float numpy arrays, one per group of
            spots, containing the desired property for each element.
    """
    all_track_property = list() # stores numpy arrays of values
    for spot_ids in all_spot_ids:
        curr_prop_values = np.empty(len(spot_ids), dtype=float)
        for i, spot_id in enumerate(spot_ids):
            spot_ind = spot_id[1]
            curr_prop_values[i] = spots[spot_ind][prop_name]    
        all_track_property.append(curr_prop_values)
    return all_track_property

def build_spot_lists_by_track(tracks_raw, spots):
    """Group spots based on the tracks they're assigned to.
    
    Args:
        tracks_raw (NodeList): DOM node list of all tracks.
        spots (list): List of dictionaries containing the attributes of each 
            spot.
            
    Returns:
        processed_tracks (tuple): Contains the following lists.
            all_track_coords (list): numpy arrays of x,y,t spot positions
            all_track_edges (list): dicts of attributes of each edge
            all_track_spot_ids (list): spot IDs and spot indices for each track
    """
    
    
    # Build x,y, t spot coordinates for each track
    # 'all_tracks' and 'spots' are lists of dictionaries directly
    # from the output of the TrackMate plugin
    
    # Build an array of spot indices
    spot_ids = np.empty([len(spots)], dtype=int)
    for i, spot in enumerate(spots):
        spot_ids[i] = int(spot['ID'])
     
    all_track_edges = list()
    all_track_coords = list() # stores numpy arrays of x,y,t coords
    all_track_spot_ids = list() # stores spot IDs and indices for each track
    for track_raw in tracks_raw:
        curr_edges = track_raw.getElementsByTagName('Edge')
        track_name = track_raw.getAttribute('name')
        edges = list()
        for edge in curr_edges:
            edge_dict = dict(edge.attributes.items())
            edge_dict['track_name'] = track_name
            edges.append(edge_dict)

        edges.sort(key=lambda k: float(k['EDGE_TIME']))
        # build spot list
        num_spots = len(edges) + 1
        track_coords = np.empty([num_spots, 3], dtype=float) # x,y,t
        track_spot_ids = np.empty([num_spots, 2], dtype=int)
        for i in range(num_spots):
            if i==(num_spots-1): # get coords of last spot
                edge = edges[i-1]
                spot_id = int(edge['SPOT_TARGET_ID'])
            else:
                edge = edges[i]
                spot_id = int(edge['SPOT_SOURCE_ID'])
            spot_ind = np.where(spot_ids==spot_id)[0][0]

            track_spot_ids[i][0] = spot_id
            track_spot_ids[i][1] = spot_ind

            track_coords[i][0] = spots[spot_ind]['POSITION_X']
            track_coords[i][1] = spots[spot_ind]['POSITION_Y']
            track_coords[i][2] = spots[spot_ind]['POSITION_T']

        all_track_edges.append(edges)
        all_track_coords.append(track_coords)
        all_track_spot_ids.append(track_spot_ids)
    
    processed_tracks = (all_track_coords, all_track_edges, 
                            all_track_spot_ids)
    return processed_tracks