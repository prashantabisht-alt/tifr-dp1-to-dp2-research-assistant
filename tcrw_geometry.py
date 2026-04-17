"""
TCRW Geometry Engine: Flexible lattice masks for arbitrary boundary conditions
===============================================================================

Supports:
  - Simple rectangular OBC (backward compatible with existing code)
  - Rectangular with internal defects (removed sites → internal boundaries)
  - Rectangular with holes (removed rectangular blocks)
  - Hybrid BC: periodic in one direction, open in the other
  - Fully custom geometries via binary array

Used by: build_transition_matrix_generic, simulate_tcrw_geometry, exact_currents_generic

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np

# Direction vectors: d=0(↑), 1(→), 2(↓), 3(←)
DX = np.array([0, 1, 0, -1], dtype=np.int64)
DY = np.array([1, 0, -1, 0], dtype=np.int64)


class LatticeMask:
    """
    Base class for lattice geometry.

    A mask defines which sites exist and how boundary conditions work.
    Subclasses override is_valid() and neighbor() for specific geometries.
    """

    def __init__(self, Lx, Ly=None):
        self.Lx = Lx
        self.Ly = Ly if Ly is not None else Lx
        self._build_index_maps()

    def is_valid(self, x, y):
        """Return True if (x, y) is a valid lattice site."""
        return 0 <= x < self.Lx and 0 <= y < self.Ly

    def neighbor(self, x, y, d):
        """
        Given site (x,y) and direction d, return (nx, ny) after move.
        Returns None if the move is blocked (OBC wall or invalid site).
        """
        nx = x + DX[d]
        ny = y + DY[d]
        if self.is_valid(nx, ny):
            return (nx, ny)
        return None

    def _build_index_maps(self):
        """Build (x,y) <-> linear index mappings for all valid sites."""
        self._sites = []
        self._site_to_idx = {}
        idx = 0
        for y in range(self.Ly):
            for x in range(self.Lx):
                if self.is_valid(x, y):
                    self._sites.append((x, y))
                    self._site_to_idx[(x, y)] = idx
                    idx += 1
        self._n_sites = idx

    @property
    def n_sites(self):
        """Number of valid lattice sites."""
        return self._n_sites

    @property
    def n_states(self):
        """Total number of states = 4 * n_sites (4 director states per site)."""
        return 4 * self._n_sites

    @property
    def valid_sites(self):
        """List of (x, y) tuples for all valid sites, in index order."""
        return list(self._sites)

    def site_to_index(self, x, y):
        """Map (x, y) -> site index. Returns None if site is invalid."""
        return self._site_to_idx.get((x, y), None)

    def index_to_site(self, idx):
        """Map site index -> (x, y)."""
        return self._sites[idx]

    def state_index(self, x, y, d):
        """Map (x, y, d) -> full state index = d * n_sites + site_index."""
        si = self.site_to_index(x, y)
        if si is None:
            return None
        return d * self._n_sites + si

    def state_to_xyd(self, state_idx):
        """Map full state index -> (x, y, d)."""
        d = state_idx // self._n_sites
        si = state_idx % self._n_sites
        x, y = self._sites[si]
        return x, y, d

    def is_boundary(self, x, y):
        """
        Check if (x, y) is a boundary site (has at least one blocked neighbor).
        """
        if not self.is_valid(x, y):
            return False
        for d in range(4):
            if self.neighbor(x, y, d) is None:
                return True
        return False

    def boundary_type(self, x, y):
        """
        Classify boundary: 'external', 'internal', 'bulk', or 'invalid'.

        External: site is on the outer rectangle boundary.
        Internal: site borders a defect/hole but not the outer boundary.
        Bulk: all 4 neighbors exist.
        """
        if not self.is_valid(x, y):
            return 'invalid'

        on_outer = (x == 0 or x == self.Lx - 1 or y == 0 or y == self.Ly - 1)
        has_blocked = any(self.neighbor(x, y, d) is None for d in range(4))

        if has_blocked and on_outer:
            return 'external'
        elif has_blocked and not on_outer:
            return 'internal'
        else:
            return 'bulk'

    def get_boundary_sites(self, btype='all'):
        """
        Return list of (x, y) boundary sites.
        btype: 'all', 'external', 'internal'
        """
        result = []
        for x, y in self._sites:
            bt = self.boundary_type(x, y)
            if btype == 'all' and bt in ('external', 'internal'):
                result.append((x, y))
            elif bt == btype:
                result.append((x, y))
        return result

    def to_array(self):
        """Return Lx × Ly boolean array. True = valid site."""
        arr = np.zeros((self.Lx, self.Ly), dtype=bool)
        for x, y in self._sites:
            arr[x, y] = True
        return arr

    def visualize(self, ax=None, show_types=True):
        """Plot the lattice mask. Color by boundary type if show_types=True."""
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(max(6, self.Lx * 0.6),
                                                    max(6, self.Ly * 0.6)))

        colors = {'bulk': '#d4e6f1', 'external': '#f5b041',
                  'internal': '#e74c3c', 'invalid': '#2c3e50'}

        # Draw grid
        for x in range(self.Lx):
            for y in range(self.Ly):
                if self.is_valid(x, y):
                    bt = self.boundary_type(x, y) if show_types else 'bulk'
                    rect = mpatches.Rectangle((x - 0.5, y - 0.5), 1, 1,
                                               facecolor=colors[bt],
                                               edgecolor='gray', lw=0.5)
                    ax.add_patch(rect)
                else:
                    rect = mpatches.Rectangle((x - 0.5, y - 0.5), 1, 1,
                                               facecolor=colors['invalid'],
                                               edgecolor='gray', lw=0.5, alpha=0.3)
                    ax.add_patch(rect)

        ax.set_xlim(-0.6, self.Lx - 0.4)
        ax.set_ylim(-0.6, self.Ly - 0.4)
        ax.set_aspect('equal')
        ax.set_xlabel('x')
        ax.set_ylabel('y')

        if show_types:
            legend_elements = [
                mpatches.Patch(facecolor=colors['bulk'], label='Bulk'),
                mpatches.Patch(facecolor=colors['external'], label='External boundary'),
                mpatches.Patch(facecolor=colors['internal'], label='Internal boundary'),
                mpatches.Patch(facecolor=colors['invalid'], alpha=0.3, label='Blocked/removed'),
            ]
            ax.legend(handles=legend_elements, loc='upper right', fontsize=8)

        return ax


class RectangleMask(LatticeMask):
    """Simple L × L rectangle with OBC on all sides. Backward compatible."""

    def __init__(self, L):
        super().__init__(L, L)


class RectangleWithDefects(LatticeMask):
    """
    L × L rectangle with some sites removed (internal defects/walls).

    Parameters
    ----------
    L : int
        Grid size.
    blocked : list of (x, y) tuples
        Sites that are removed from the lattice.
    """

    def __init__(self, L, blocked):
        self._blocked = set(blocked)
        super().__init__(L, L)

    def is_valid(self, x, y):
        if not (0 <= x < self.Lx and 0 <= y < self.Ly):
            return False
        return (x, y) not in self._blocked


class RectangleWithHoles(LatticeMask):
    """
    L × L rectangle with rectangular holes removed.

    Parameters
    ----------
    L : int
        Grid size.
    holes : list of (x0, y0, w, h) tuples
        Each hole removes sites in [x0, x0+w) × [y0, y0+h).
    """

    def __init__(self, L, holes):
        self._blocked = set()
        for x0, y0, w, h in holes:
            for x in range(x0, x0 + w):
                for y in range(y0, y0 + h):
                    self._blocked.add((x, y))
        super().__init__(L, L)

    def is_valid(self, x, y):
        if not (0 <= x < self.Lx and 0 <= y < self.Ly):
            return False
        return (x, y) not in self._blocked


class HybridBCMask(LatticeMask):
    """
    L × L lattice with per-axis boundary conditions.

    Parameters
    ----------
    L : int
        Grid size.
    pbc_x : bool
        If True, x-direction is periodic.
    pbc_y : bool
        If True, y-direction is periodic.
    """

    def __init__(self, L, pbc_x=False, pbc_y=False):
        self.pbc_x = pbc_x
        self.pbc_y = pbc_y
        super().__init__(L, L)

    def neighbor(self, x, y, d):
        nx = x + DX[d]
        ny = y + DY[d]

        # Handle PBC wrapping
        if self.pbc_x:
            nx = nx % self.Lx
        if self.pbc_y:
            ny = ny % self.Ly

        if self.is_valid(nx, ny):
            return (nx, ny)
        return None

    def boundary_type(self, x, y):
        """With PBC axes, only OBC axes contribute to boundary."""
        if not self.is_valid(x, y):
            return 'invalid'
        has_blocked = any(self.neighbor(x, y, d) is None for d in range(4))
        if has_blocked:
            return 'external'
        return 'bulk'


class PBCWithDefects(LatticeMask):
    """
    L × L lattice with full PBC on outer boundary + internal defects.

    This is for Fig 11(c)-(d): outer boundary wraps (PBC), but there
    are removed sites in the interior creating internal boundaries.
    """

    def __init__(self, L, blocked):
        self._blocked = set(blocked)
        super().__init__(L, L)

    def is_valid(self, x, y):
        # Any position in [0, L) is geometrically valid...
        if not (0 <= x < self.Lx and 0 <= y < self.Ly):
            return False
        # ...unless it's a blocked/defect site
        return (x, y) not in self._blocked

    def neighbor(self, x, y, d):
        # PBC wrapping on outer boundary
        nx = (x + DX[d]) % self.Lx
        ny = (y + DY[d]) % self.Ly
        if self.is_valid(nx, ny):
            return (nx, ny)
        return None  # blocked by defect

    def boundary_type(self, x, y):
        if not self.is_valid(x, y):
            return 'invalid'
        has_blocked = any(self.neighbor(x, y, d) is None for d in range(4))
        if has_blocked:
            return 'internal'  # PBC has no external boundary
        return 'bulk'


class CustomMask(LatticeMask):
    """
    Build from a 2D boolean array. True = valid site.
    Optionally specify per-axis PBC.
    """

    def __init__(self, grid, pbc_x=False, pbc_y=False):
        self._grid = np.asarray(grid, dtype=bool)
        self.pbc_x = pbc_x
        self.pbc_y = pbc_y
        Lx, Ly = self._grid.shape
        super().__init__(Lx, Ly)

    def is_valid(self, x, y):
        xx = x % self.Lx if self.pbc_x else x
        yy = y % self.Ly if self.pbc_y else y
        if not (0 <= xx < self.Lx and 0 <= yy < self.Ly):
            return False
        return bool(self._grid[xx, yy])

    def neighbor(self, x, y, d):
        nx = x + DX[d]
        ny = y + DY[d]
        if self.pbc_x:
            nx = nx % self.Lx
        if self.pbc_y:
            ny = ny % self.Ly
        if self.is_valid(nx, ny):
            return (nx, ny)
        return None


# ============================================================
# Predefined defect geometries for the paper's figures
# ============================================================

def paper_defect_edge_notch(L=10):
    """
    Fig 2(i)-(k): Edge deformation — notch cut into outer boundary.
    Remove a rectangular notch from the right edge.
    """
    blocked = []
    # Remove 3×3 block from right edge, centered vertically
    cx, cy = L - 1, L // 2
    for dx in range(3):
        for dy in range(-1, 2):
            bx = cx - dx
            by = cy + dy
            if 0 <= bx < L and 0 <= by < L:
                blocked.append((bx, by))
    return RectangleWithDefects(L, blocked)


def paper_defect_internal_block(L=10):
    """
    Fig 2(l)-(o): Internal defect — removed block in the bulk.
    Creates an internal boundary with its own edge current.
    """
    # 3×3 hole in the center
    cx, cy = L // 2, L // 2
    blocked = []
    for dx in range(-1, 2):
        for dy in range(-1, 2):
            blocked.append((cx + dx, cy + dy))
    return RectangleWithDefects(L, blocked)


def paper_defect_two_holes(L=15):
    """
    Fig 12: Two separated internal defects creating disconnected internal boundaries.
    """
    blocked = []
    # Hole 1: 3×3 at (4, 4)
    for dx in range(3):
        for dy in range(3):
            blocked.append((3 + dx, 3 + dy))
    # Hole 2: 3×3 at (10, 10)
    for dx in range(3):
        for dy in range(3):
            blocked.append((10 + dx, 10 + dy))
    return RectangleWithDefects(L, blocked)


def paper_pbc_with_internal_defect(L=10):
    """
    Fig 11(c)-(d): PBC outer boundary with internal defect block.
    """
    cx, cy = L // 2, L // 2
    blocked = []
    for dx in range(-1, 2):
        for dy in range(-1, 2):
            blocked.append((cx + dx, cy + dy))
    return PBCWithDefects(L, blocked)


# ============================================================
# Quick self-test
# ============================================================

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    masks = [
        ("Rectangle L=8", RectangleMask(8)),
        ("Edge notch L=10", paper_defect_edge_notch(10)),
        ("Internal block L=10", paper_defect_internal_block(10)),
        ("Two holes L=15", paper_defect_two_holes(15)),
        ("Hybrid (PBC-y, OBC-x) L=8", HybridBCMask(8, pbc_x=False, pbc_y=True)),
        ("PBC + internal defect L=10", paper_pbc_with_internal_defect(10)),
    ]

    for ax, (name, mask) in zip(axes.flat, masks):
        mask.visualize(ax=ax)
        ax.set_title(f"{name}\n({mask.n_sites} sites, {mask.n_states} states)", fontsize=10)

    plt.tight_layout()
    plt.savefig('tcrw_geometry_test.png', dpi=150)
    plt.close()
    print("Saved tcrw_geometry_test.png")

    # Verify backward compatibility
    rect = RectangleMask(10)
    assert rect.n_sites == 100
    assert rect.n_states == 400
    assert rect.state_index(3, 5, 2) == 2 * 100 + 5 * 10 + 3
    print("RectangleMask backward compatibility: OK")

    # Verify defect mask
    defect = paper_defect_internal_block(10)
    assert defect.n_sites == 100 - 9  # 9 sites removed
    assert defect.is_valid(5, 5) == False  # center blocked
    assert defect.is_valid(0, 0) == True
    int_bnd = defect.get_boundary_sites('internal')
    ext_bnd = defect.get_boundary_sites('external')
    print(f"Internal block: {len(int_bnd)} internal boundary sites, {len(ext_bnd)} external")

    # Verify PBC with defect
    pbc_def = paper_pbc_with_internal_defect(10)
    assert pbc_def.boundary_type(0, 0) == 'bulk'  # outer edge is PBC, not boundary
    assert pbc_def.boundary_type(3, 5) == 'internal'  # next to defect
    print(f"PBC+defect: {len(pbc_def.get_boundary_sites('internal'))} internal boundary sites")

    print("\nAll geometry tests passed.")
