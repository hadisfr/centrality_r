#!/usr/bin/env Rscript

# Algorithmic Graph Theory ~ F 97 :: CA2 ~ Q2
# by Hadi Safari (hadi@hadisafari.ir)

library(igraph)
library(R2HTML)
gnet <- read_graph("gnet.graphml", format = "graphml")

# We could use `decompose` to make connected components seperate

################

make_nodes_table <- function(g) {

    get_top_10_nodes <- function(c_by_node_index) {
        return(paste(V(g)[sort(c_by_node_index, decreasing = TRUE, index.return = TRUE)$ix[c(1:10)]]$name, collapse = ", "))
    }

    n <- length(V(g))
    u <- as.undirected(g)

    df <- data.frame()

    methods <- c("alpha_centrality", "authority_score", "betweenness", "closeness",
        # "rwalk_closeness",
        "eigen_centrality", "hub_score", "subgraph_centrality", "power_centrality")

    df <- rbind(
        c(get_top_10_nodes(alpha_centrality(g)), get_top_10_nodes(alpha_centrality(u))),
        c(get_top_10_nodes(authority_score(g)$vector), get_top_10_nodes(authority_score(u)$vector)),
        c(get_top_10_nodes(betweenness(g)), get_top_10_nodes(betweenness(u, directed = FALSE))),
        c(get_top_10_nodes(closeness(g)), get_top_10_nodes(closeness(u))),
        # c(get_top_10_nodes(rwalk_closeness(g)), get_top_10_nodes(rwalk_closeness(u))),
        c(get_top_10_nodes(eigen_centrality(g)$vector), get_top_10_nodes(eigen_centrality(u, directed = FALSE)$vector)),
        c(get_top_10_nodes(hub_score(g)$vector), get_top_10_nodes(hub_score(u)$vector)),
        c(get_top_10_nodes(subgraph_centrality(g)), get_top_10_nodes(subgraph_centrality(u))),
        c(get_top_10_nodes(power_centrality(g)), get_top_10_nodes(power_centrality(u)))
    )

    colnames(df) <- c("Top 10 nodes in directed graph", "Top 10 nodes in undirected graph")
    rownames(df) <- methods

    html = HTMLInitFile(getwd(), filename = "top10", CSSFile = "top10table.css")
    HTML(as.title("Top 10 Nodes of 'gnet' by Closeness"), file = html)
    HTML(df, file = html)
    browseURL(paste("file://", html))

    return(df)
}

########

rwalk_closeness <- function(g) {
    # based on implementation of tidygraph and netrankr packages
    n <- vcount(g)
    A <- get.adjacency(g, sparse = FALSE)
    M <- A / rowSums(A)
    e <- rep(1, n - 1)
    H <- matrix(0, n, n)
    lapply(1:n, function(j) H[j, -j] <<- solve(diag(e) - M[-j, -j]) %*% e)
    # Or:
        # for (j in 1:n) {
        #     Mj <- M[-j, -j]
        #     Hij <- solve(diag(e) - Mj) %*% e
        #     H[j, -j] <- Hij # transposed to original to fit framework with rowSums
        # }
    if(!is.null(V(g)$name))
        rownames(H) <- colnames(H) <- V(g)$name
    return(rowSums(H)^-1)
}

################

make_nodes_table(gnet)

########

c <- closeness(gnet)
rwc <- rwalk_closeness(gnet)

print("closeness:")
print(c)
print("random walk closeness:")
print(rwc)
