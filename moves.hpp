/*
 * - MK : add/remove empty cluster
 *
 * Partitioning of edges:
 * - Gibbs update of one edge
 * - SM-like move
 * - M3-like move
 *
 * Overlapping edges.
 * - Deprecated: M-H step on one edge.  Propose for each community independently, accept if and only if the row is non-empty
 * - Surely there is a Gibbs update I can do?
 */