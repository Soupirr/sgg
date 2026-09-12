"""Convert a Biopython tree into plotly coordinates / a Figure."""

import plotly.graph_objects as go

from sgg.analyzer import get_color
from sgg.palette import PALETTE


def tree_to_plotly(tree, mode: str = "cladogram"):
    x_lines, y_lines, x_nodes, y_nodes, labels = [], [], [], [], []
    y_pos = {}
    counter = [0]

    def get_y(clade):
        if clade.is_terminal():
            y_pos[clade] = counter[0]
            counter[0] += 1
        else:
            for c in clade.clades:
                get_y(c)
            y_pos[clade] = sum(y_pos[c] for c in clade.clades) / len(clade.clades)

    def get_x(clade, x=0):
        clade.x = x
        for c in clade.clades:
            if mode == "phylogram":
                get_x(c, x + (c.branch_length or 0))
            else:
                get_x(c, x + 1)

    get_y(tree.root)
    get_x(tree.root)

    max_x = max(c.x for c in tree.find_clades())

    def collect(clade):
        for c in clade.clades:
            x_lines.extend([clade.x, c.x, None])
            y_lines.extend([y_pos[clade], y_pos[clade], None])
            x_lines.extend([c.x, c.x, None])
            y_lines.extend([y_pos[clade], y_pos[c], None])
            collect(c)

        if clade.is_terminal():
            x_lines.extend([clade.x, clade.x + 0.5, None])
            y_lines.extend([y_pos[clade], y_pos[clade], None])
            x_nodes.append(clade.x + 0.5)
            labels.append((clade.name or "")[:70])
        else:
            x_nodes.append(clade.x)
            confidence = clade.confidence
            if confidence is not None and confidence >= 0.5:
                labels.append(f"{confidence * 100:.1f}")
            else:
                labels.append("")

        y_nodes.append(y_pos[clade])

    collect(tree.root)
    return x_lines, y_lines, x_nodes, y_nodes, labels, counter[0], max_x


def build_tree_figure(tree, title: str, mode: str = "cladogram", multi_query: bool = False):
    x_lines, y_lines, x_nodes, y_nodes, labels, n_leaves, _max_x = tree_to_plotly(tree, mode=mode)

    if multi_query:
        q_labels = list(dict.fromkeys(lab for lab in labels if lab.startswith("QUERY_")))
        q_colors = {lab: PALETTE[i % len(PALETTE)] for i, lab in enumerate(q_labels)}
    else:
        q_colors = {}

    node_colors = [
        q_colors.get(lab, "#00FF00") if lab.startswith("QUERY_") else get_color(lab)
        for lab in labels
    ]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x_lines,
            y=y_lines,
            mode="lines",
            line={"color": "rgba(150,150,150,0.3)", "width": 1.5},
            hoverinfo="none",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x_nodes,
            y=y_nodes,
            mode="markers+text",
            marker={"size": 6, "color": node_colors},
            text=labels,
            textposition="middle right",
            textfont={
                "size": 12,
                "color": [
                    q_colors.get(lab, "#2ECC71") if lab.startswith("QUERY_") else "white"
                    for lab in labels
                ],
            },
            hoverinfo="text",
        )
    )
    fig.update_layout(
        title={"text": title[:66], "font": {"size": 24, "color": "white"}, "x": 0, "xanchor": "left"},
        showlegend=False,
        height=max(400, n_leaves * 25),
        width=1800,
        margin={"l": 0, "r": 0, "t": 40, "b": 0},
        xaxis={"showgrid": False, "zeroline": False, "showticklabels": False, "fixedrange": False},
        yaxis={"showgrid": False, "zeroline": False, "showticklabels": False, "fixedrange": False},
        plot_bgcolor="#060d14",
        paper_bgcolor="rgba(0,0,0,0)",
        shapes=[
            {
                "type": "rect",
                "xref": "paper",
                "yref": "paper",
                "x0": 0,
                "y0": 0,
                "x1": 1,
                "y1": 1,
                "line": {"color": "rgba(255,255,255,0.3)", "width": 1},
            }
        ],
        dragmode="pan",
    )
    return fig
