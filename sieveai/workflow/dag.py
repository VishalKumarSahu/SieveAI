"""Workflow DAG - dependency resolution and topological sort."""

from collections import defaultdict, deque
from typing import Dict, List, Any


class WorkflowDAG:
    """
    Build execution order from plugin dependencies.
    
    Plugins declare 'requires' and 'provides'.
    DAG auto-sorts to satisfy dependencies

    Example:
        structuresync: requires=[], provides=[receptors, ligands]
        vina: requires=[structuresync], provides=[docked_complexes]
        chimerax: requires=[vina], provides=[analysis]
        
        Order: structuresync → vina → chimerax
    """
    
    def __init__(self, plugins: Dict[str, Dict[str, Any]]):
        self.plugins = plugins
        self.graph = defaultdict(list)
        self.indegree = defaultdict(int)
    
    def build(self) -> List[str]:
        """
        Build dependency graph, return topological order.
        
        Returns:
            List of plugin uids in execution order.
        
        Raises:
            ValueError: If circular dependency detected.
        """
        # Build graph from 'requires' declarations
        for uid, info in self.plugins.items():
            requires = info.get("requires", [])
            for req in requires:
                self.graph[req].append(uid)
                self.indegree[uid] += 1
            
            # Initialize indegree for plugins with no requirements
            if uid not in self.indegree:
                self.indegree[uid] = 0
        
        # Topological sort (Kahn's algorithm)
        queue = deque([uid for uid in self.plugins if self.indegree[uid] == 0])
        order = []
        
        while queue:
            uid = queue.popleft()
            order.append(uid)
            
            for dependent in self.graph[uid]:
                self.indegree[dependent] -= 1
                if self.indegree[dependent] == 0:
                    queue.append(dependent)
        
        # Detect cycles
        if len(order) != len(self.plugins):
            remaining = [uid for uid in self.plugins if uid not in order]
            raise ValueError(
                f"Circular dependency detected in plugins: {remaining}"
            )
        
        return order
    
    def get_dependencies(self, uid: str) -> List[str]:
        """Get all plugins that must run before uid."""
        deps = []
        for other_uid, info in self.plugins.items():
            if uid in info.get("requires", []):
                deps.append(other_uid)
        return deps
    
    def get_dependents(self, uid: str) -> List[str]:
        """Get all plugins that depend on uid."""
        return self.graph.get(uid, [])
