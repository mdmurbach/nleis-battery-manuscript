import numpy as np
from cycler import cycler
from matplotlib.ticker import ScalarFormatter


class FixedOrderFormatter(ScalarFormatter):
    """Formats axis ticks using scientific notation with a constant order of
    magnitude"""
    def __init__(self, order_of_mag=0, useOffset=True, useMathText=True):
        self._order_of_mag = order_of_mag
        ScalarFormatter.__init__(self, useOffset=useOffset,
                                 useMathText=useMathText)

    def _set_orderOfMagnitude(self, range):
        """Over-riding this to avoid having orderOfMagnitude reset elsewhere"""
        self.orderOfMagnitude = self._order_of_mag


def plot_nyquist(fig, ax, frequencies, Zs, harmonic, subplot_label,
                 limits, mask=None, names=[], markers=[]):

    ax.set_prop_cycle(cycler('color', ['#1f77b4', '#ff7f0e',
                                       '#2ca02c', '#d62728', '#9467bd']) +
                      cycler('ls', ['dashed', 'solid',
                                    'dashdot', 'dotted', (0, (6, 2))]))

    for i, Z in enumerate(Zs):

        props = next(ax._get_lines.prop_cycler)
        color = props['color']
        ls = props['linestyle']

        if mask is not None:
            Z = Z[mask]

        if len(names) > 0:
            ax.plot(np.real(Z), -np.imag(Z), linewidth=2.5, color=color, ls=ls,
                    label=names[i])
        else:
            ax.plot(np.real(Z), -np.imag(Z), linewidth=2.5, color=color, ls=ls)

        mk_symbols = 's^vs^v'
        if len(markers) > 1:
            for j, mk_freq in enumerate(markers):
                if mk_freq is not None:
                    mk_loc = np.where(np.array(frequencies) == mk_freq)[0][0]
                    ax.plot(np.real(Z[mk_loc]), -np.imag(Z[mk_loc]),
                            mk_symbols[j], markersize=8, color=color,
                            mec='k', zorder=10)

    if harmonic == 1:
        ax.set_xlabel(r'$\tilde{Z}_{1,1}^{\prime}(\omega)$ [$\Omega$]',
                      fontsize=18)
        ax.set_ylabel(r'$-\tilde{Z}_{1,1}^{\prime\prime}(\omega)$ [$\Omega$]',
                      fontsize=18)
    elif harmonic == 2:
        ax.set_xlabel(r'$\tilde{Z}_{2,2}^{\prime}(\omega)$ [$\Omega / A$]',
                      fontsize=18)
        ax.set_ylabel(r'$-\tilde{Z}_{2,2}^{\prime\prime}(\omega)$' +
                      r' [$\Omega / A$]', fontsize=18)

    ax.set_xlim(limits[0])
    ax.set_ylim(limits[1])
    ax.xaxis.set_major_formatter(FixedOrderFormatter(limits[2]))
    ax.yaxis.set_major_formatter(FixedOrderFormatter(limits[2]))
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.locator_params(nbins=5)
    ax.grid(b=True, which='major', axis='both', alpha=.5)
    y_offset = ax.yaxis.get_offset_text()
    y_offset.set_size(16)

    t = ax.xaxis.get_offset_text()
    t.set_size(16)
    t.set_position((1.08, 10))

    ax.text(.95, .95, subplot_label, weight='bold', family='sans-serif',
            fontsize=28, horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes)

    if len(names) > 1:
        ax.legend(loc='upper left', frameon=True, fontsize=16)

    return ax


def plot_nyquist_sim(fig, ax, frequencies, Zs, harmonic, subplot_label,
                     limits, mask=None, names=[], markers=[]):

    ax.set_prop_cycle(cycler('color', ['#1f77b4', '#ff7f0e',
                                       '#2ca02c', '#d62728', '#9467bd']) +
                      cycler('ls', ['dashed', 'solid',
                                    'dashdot', 'dotted', (0, (6, 2))]))

    for i, Z in enumerate(Zs):

        props = next(ax._get_lines.prop_cycler)
        color = props['color']
        ls = props['linestyle']

        if mask is not None:
            Z = Z[mask]

        if len(names) > 0:
            ax.plot(np.real(Z), -np.imag(Z), linewidth=2.5, color=color, ls=ls,
                    label=names[i])
        else:
            ax.plot(np.real(Z), -np.imag(Z), linewidth=2.5, color=color, ls=ls)

        mk_symbols = 's^vs^v'
        if len(markers) > 1:
            for j, mk_freq in enumerate(markers):
                if mk_freq is not None:
                    mk_loc = np.where(np.array(frequencies) == mk_freq)[0][0]
                    ax.plot(np.real(Z[mk_loc]), -np.imag(Z[mk_loc]),
                            mk_symbols[j], markersize=8, color=color,
                            mec='k', zorder=10)

    if harmonic == 1:
        ax.set_xlabel(r'$\tilde{Z}_{1,1}^{\prime}(\omega)$ [$\Omega$ - $m^2$]',
                      fontsize=16)
        ax.set_ylabel(r'$-\tilde{Z}_{1,1}^{\prime\prime}(\omega)$' +
                      ' [$\Omega$ - $m^2$]', fontsize=16)
    elif harmonic == 2:
        ax.set_xlabel(r'$\tilde{Z}_{2,2}^{\prime}(\omega)$' +
                      ' [$\Omega / A$ - $m^4$]', fontsize=16)
        ax.set_ylabel(r'$-\tilde{Z}_{2,2}^{\prime\prime}(\omega)$' +
                      ' [$\Omega / A$ - $m^4$]', fontsize=16)

        label = ax.xaxis.get_label()
        x_lab_pos, y_lab_pos = label.get_position()
        label.set_position([0.1, y_lab_pos])
        label.set_horizontalalignment('left')
        ax.xaxis.set_label(label)

    ax.set_xlim(limits[0])
    ax.set_ylim(limits[1])
    ax.xaxis.set_major_formatter(FixedOrderFormatter(limits[2]))
    ax.yaxis.set_major_formatter(FixedOrderFormatter(limits[2]))
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.locator_params(nbins=5)
    ax.grid(b=True, which='major', axis='both', alpha=.5)
    y_offset = ax.yaxis.get_offset_text()
    y_offset.set_size(16)

    t = ax.xaxis.get_offset_text()
    t.set_size(16)
    t.set_position((1.08, 10))

    ax.text(.95, .95, subplot_label, weight='bold', family='sans-serif',
            fontsize=28, horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes)

    if len(names) > 1:
        ax.legend(loc='upper left', frameon=True, fontsize=16)

    return ax
