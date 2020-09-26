
import re
import os, sys

import numpy as np
import pandas as pd
import networkx as nx

from networkx.drawing.nx_agraph import graphviz_layout


def add_graphics_theme_to_nx_graph(
        nx_graph,
        edge_color=None,
        node_size_factor=50,
        edge_size_factor=500):
    """adjust nodes and edges
    """
    # node size, stroke
    for node_name, node_attrs in nx_graph.nodes(data=True):

        #node_size = nx_graph.nodes[node_name]["numexamples"] / float(node_size_factor)
        #node_size = nx_graph.nodes[node_name]["numexamples"] / float(nx_graph.graph["numexamples"])
        #node_size *= node_size_factor
        node_size = 100
        
        graphics = {
            "type": "ellipse",
            "w": node_size,
            "h": node_size,
            "fill": "#FFFFFF",
            "outline": "#000000",
            "width": 1.0,
            "fontSize": 14
        }

        if nx_graph.nodes[node_name].get("graphics") is not None:
            nx_graph.nodes[node_name]["graphics"].update(graphics)
        else:
            nx_graph.nodes[node_name]["graphics"] = graphics

    # edges
    for start_node, end_node in nx_graph.edges():
        for edge_idx in xrange(len(nx_graph[start_node][end_node])):

            #edge_width = nx_graph[start_node][end_node][edge_idx]["numexamples"] / float(
            #    edge_size_factor)
            
            #edge_width = nx_graph[start_node][end_node][edge_idx]["numexamples"] / float(
            #    nx_graph.graph["numexamples"])
            #edge_width *= edge_size_factor 
            edge_width = 1.0
            
            graphics = {
                "type": "arc",
                "width": edge_width,
                "targetArrow": "delta"
            }
            
            if edge_color is not None:
                graphics["fill"] = edge_color

            if nx_graph[start_node][end_node][edge_idx].get("graphics") is not None:
                nx_graph[start_node][end_node][edge_idx]["graphics"].update(graphics)
            else:
                nx_graph[start_node][end_node][edge_idx]["graphics"] = graphics
            
    return None


def stringize_nx_graph(nx_graph):
    """preparatory function for writing out to gml
    """
    # graph attributes
    for key in nx_graph.graph.keys():
        if isinstance(nx_graph.graph[key], (list, set, np.ndarray)):
            nx_graph.graph[key] = ",".join([
                str(val) for val in list(nx_graph.graph[key])])

    # node attributes
    for node_name, node_attrs in nx_graph.nodes(data=True):
        for key in node_attrs.keys():
            if isinstance(nx_graph.nodes[node_name][key], (list, set, np.ndarray)):
                nx_graph.nodes[node_name][key] = ",".join([
                    str(val) for val in nx_graph.nodes[node_name][key]])
        # adjust node name for nice output in cytoscape
        new_node_name = re.sub(r"HCLUST.\d+_", "", node_name)
        new_node_name = new_node_name.replace(".UNK.0.A", "")
        nx_graph.nodes[node_name]["name"] = new_node_name
                
    # edge attributes
    for start_node, end_node in nx_graph.edges():
        for edge_idx in xrange(len(nx_graph[start_node][end_node])):
            edge_attrs = nx_graph[start_node][end_node][edge_idx]
            for key in edge_attrs.keys():
                if isinstance(edge_attrs[key], (list, set, np.ndarray)):
                    nx_graph[start_node][end_node][edge_idx][key] = ",".join([
                        str(val) for val in nx_graph[start_node][end_node][edge_idx][key]])
                    
    return nx_graph


def main():
    """build a network view
    """
    # files
    summary_file = sys.argv[1]
    pwms_to_tfs_file = sys.argv[2]
    expressed_tfs_file = sys.argv[3] # TODO

    # TODO pull in num regions to resize things? but complicated with overlaps etc
    # TODO edit edges with type of interaction
    # TODO may want to color by trajectory, to demonstrate waves of trajectory
    
    # read in data
    summary = pd.read_csv(summary_file, sep="\t")
    pwms_to_tfs = pd.read_csv(pwms_to_tfs_file, sep="\t")
    pwms_to_tfs = pwms_to_tfs[pwms_to_tfs["expressed"].notna()]
    pwms_to_filt_tfs = {} # dict: key - pwm names, vals - dict of ensembl ids to hgnc ids
    for line_idx in range(pwms_to_tfs.shape[0]):
        pwm_info = pwms_to_tfs.iloc[line_idx,:]
        pwm_name = pwm_info["hclust_model_name"]
        pwm_to_tf = dict(zip(pwm_info["expressed"].split(";"), pwm_info["expressed_hgnc"].split(";")))
        pwms_to_filt_tfs[pwm_name] = pwm_to_tf

        
    # filter expressed hgncs for dynamic ones only
    tfs_filt = pd.read_csv(expressed_tfs_file, sep="\t", index_col=0)
    for pwm_name in pwms_to_filt_tfs.keys():
        tfs_tmp = pwms_to_filt_tfs[pwm_name]
        for ensembl_tf in tfs_tmp.keys():
            if ensembl_tf not in tfs_filt.index:
                del tfs_tmp[ensembl_tf]
        if len(tfs_tmp.keys()) == 0:
            del pwms_to_filt_tfs[pwm_name]
        pwms_to_filt_tfs[pwm_name] = tfs_tmp

    # add in tfs column
    tf1 = []
    for pwm in summary["pwm1"]:
        tf_str = []
        for ensembl_id in pwms_to_filt_tfs[pwm]:
            tf_str.append(pwms_to_filt_tfs[pwm][ensembl_id])
            # TODO try add in max point
            expression = tfs_filt.loc[ensembl_id,:]
            max_idx = np.argmax(expression.values)
            tf_str.append(str(max_idx))
        tf_str = (";").join(tf_str)
        tf1.append(tf_str)
    summary["tf1"] = tf1

    tf2 = []
    for pwm in summary["pwm2"]:
        tf_str = []
        for ensembl_id in pwms_to_filt_tfs[pwm]:
            tf_str.append(pwms_to_filt_tfs[pwm][ensembl_id])
            expression = tfs_filt.loc[ensembl_id,:]
            max_idx = np.argmax(expression.values)
            tf_str.append(str(max_idx))
        tf_str = (";").join(tf_str)
        tf2.append(tf_str)
    summary["tf2"] = tf2
    
    # remove failed rules
    summary = summary[~summary["interaction"].str.contains("FAILED")]
    
    # make graph
    graph = nx.from_pandas_edgelist(summary, "tf1", "tf2")

    # set up positions
    #pos = graphviz_layout(graph, prog="dot")
    pos = graphviz_layout(graph, prog="neato")
    scale_factor = 3
    for key in pos.keys():
        coords = pos[key]
        pos[key] = {"x": scale_factor*coords[0], "y": -scale_factor*coords[1]}
    nx.set_node_attributes(graph, pos, "graphics") # note this is diff from v1 to v2 in networkx
    
    # add graphics
    add_graphics_theme_to_nx_graph(graph)

    # write gml
    out_file = "summary.gml"
    nx.write_gml(stringize_nx_graph(graph), out_file, stringizer=str)

    # tfs: for each tf, get gene column
    
    
    return

main()
