# Function to generate a Sankey plot for chromatin state transitions
plot_sankey <- function(from_cell_bed, to_cell_bed, fromstate) {
  
  # Create a temporary file to store the output
  tempfile_for_output <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = "")
  print(tempfile_for_output)
  print(fromstate)
  # System command to filter data for specified state transitions and save to temporary file
  print(paste0(
    'tail -n +2 ', from_cell_bed[['datapath']],
    ' | cut -f 1,2,3,4 | ggrep -P "\\t', fromstate, '$"',
    ' | bedtools intersect -a stdin -b ', to_cell_bed[['datapath']], 
    ' -wao | ggrep -v -P "\\t\\.\\t\\d+$" | ',
    'awk \'BEGIN{OFS="\\t"} {a[$8]+=$14} END {for (i in a) print "',
    fromstate, '",i,a[i]}\'| sort -n -k2,2 > ', tempfile_for_output
  ))
  system(paste0(
    'tail -n +2 ', from_cell_bed[['datapath']],
    ' | cut -f 1,2,3,4 | ggrep -P "\\t', fromstate, '$"',
    ' | bedtools intersect -a stdin -b ', to_cell_bed[['datapath']], 
    ' -wao | ggrep -v -P "\\t\\.\\t\\d+$" | ',
    'awk \'BEGIN{OFS="\\t"} {a[$8]+=$14} END {for (i in a) print "',
    fromstate, '",i,a[i]}\'| sort -n -k2,2 > ', tempfile_for_output
  ))
  
  # Read in the processed data from the temporary file
  state_map_base_counts <- readr::read_tsv(
    tempfile_for_output,
    col_names = c("Base_State", "New_State", "Count"),
    col_types = "cci"
  )
  
  # Number of rows in the data for creating Sankey plot links
  num_row_for_sankey <- nrow(state_map_base_counts)
  
  # Define colors for plot nodes
  plot_cols <- c(
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", 
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", 
    "#A0CBE8", "#FFBE7D", "#8CD17D", "#F1CE63", "#86BCB6", 
    "#FABFD2", "#FF9D9A", "#D4A6C8", "#D7B5A6", "#BAB0AC"
  )
  
  # Display labels for the 'from' and 'to' states
  query_state_display_string <- from_cell_bed[['name']]
  db_state_display_string <- to_cell_bed[['name']]
  
  # Define sources and targets for Sankey plot links
  source_for_sankey <- rep(0, num_row_for_sankey)  # From state base
  target_for_sankey <- seq(1, num_row_for_sankey)  # To state target
  
  print(unique(state_map_base_counts$Base_State))
  # Create the Sankey plot using plotly
  fig <- plot_ly(
    type = "sankey",
    orientation = "h",
    node = list(
      label = c(
        unique(state_map_base_counts$Base_State),
        unique(state_map_base_counts$New_State)
      ),
      customdata = c(
        unique(state_map_base_counts$Base_State),
        unique(state_map_base_counts$New_State)
      ),
      hovertemplate = "State: %{label}",
      color = plot_cols,
      pad = 25,
      line = list(color = "black", width = 2)
    ),
    link = list(
      source = source_for_sankey,
      target = target_for_sankey,
      value = state_map_base_counts$Count,
      color = adjustcolor(plot_cols[2:length(plot_cols)], alpha.f = 0.2),
      hovertemplate = paste0(
        query_state_display_string, ": %{source.customdata}<br />",
        db_state_display_string, ": %{target.customdata}"
      )
    ),
    textfont = list(size = 16)
  ) %>%
    config(displaylogo = FALSE) %>%
    layout(
      title = list(text = "", automargin = TRUE),
      modebar = list(remove = c("hoverCompareCartesian")),
      autosize = TRUE,
      width = 680,
      height = 500,
      transition = list(easing = "bounce-out")
    ) %>%
    # Add annotations for state labels
    add_annotations(
      x = 1, y = 0.5, xshift = 45,
      text = db_state_display_string,
      font = list(size = 16),
      showarrow = FALSE,
      textangle = 90
    ) %>%
    add_annotations(
      x = 0, y = 0.5, xshift = -45,
      text = query_state_display_string,
      font = list(size = 16),
      showarrow = FALSE,
      textangle = -90
    )
  
  # Return the final Sankey plot object
  return(fig)
}