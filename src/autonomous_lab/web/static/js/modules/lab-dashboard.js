/**
 * Autonomous Lab Dashboard Module
 *
 * Fetches lab state from /api/autolab/state and renders
 * a dashboard panel showing role indicator, figure gallery,
 * paper progress, and file counts.
 *
 * Activates automatically when the summary contains [AUTOLAB].
 */

class LabDashboard {
    constructor() {
        this.container = null;
        this.state = null;
        this.active = false;
    }

    /**
     * Check if the current session is an Autonomous Lab session
     * by looking for the [AUTOLAB] marker in the summary.
     */
    checkAutolab(summaryText) {
        return summaryText && summaryText.includes('[AUTOLAB]');
    }

    /**
     * Initialize the dashboard. Call after page load.
     */
    async init() {
        // Check if dashboard container exists; if not, create it
        this.container = document.getElementById('lab-dashboard');
        if (!this.container) {
            this.container = document.createElement('div');
            this.container.id = 'lab-dashboard';
            this.container.className = 'lab-dashboard';
            // Insert before the feedback form
            const summarySection = document.getElementById('combinedSummaryContent') ||
                                   document.querySelector('.summary-section');
            if (summarySection && summarySection.parentNode) {
                summarySection.parentNode.insertBefore(this.container, summarySection.nextSibling);
            }
        }

        await this.refresh();
    }

    /**
     * Fetch state from the API and re-render.
     */
    async refresh() {
        try {
            const response = await fetch('/api/autolab/state');
            this.state = await response.json();
            this.active = this.state && this.state.active;

            if (this.active) {
                this.render();
                if (this.container) {
                    this.container.style.display = 'block';
                }
            } else {
                if (this.container) {
                    this.container.style.display = 'none';
                }
            }
        } catch (err) {
            console.warn('Lab dashboard fetch error:', err);
            if (this.container) {
                this.container.style.display = 'none';
            }
        }
    }

    /**
     * Render the dashboard HTML.
     */
    render() {
        if (!this.container || !this.state) return;

        const s = this.state;
        const roleClass = s.next_role === 'pi' ? 'role-pi' : 'role-trainee';
        const roleLabel = s.next_role === 'pi' ? 'PI Turn' : 'Trainee Turn';

        // Paper progress bars
        let paperHtml = '';
        if (s.paper_progress) {
            for (const [section, info] of Object.entries(s.paper_progress)) {
                const words = info.words || 0;
                const barWidth = Math.min(words / 5, 100); // rough scale: 500 words = 100%
                const statusClass = words > 0 ? 'has-content' : 'empty';
                paperHtml += `
                    <div class="paper-section ${statusClass}">
                        <span class="section-name">${section}</span>
                        <div class="progress-bar">
                            <div class="progress-fill" style="width: ${barWidth}%"></div>
                        </div>
                        <span class="word-count">${words > 0 ? words + 'w' : '--'}</span>
                    </div>
                `;
            }
        }

        // File counts
        const fc = s.file_counts || {};
        let filesHtml = '';
        for (const [dir, count] of Object.entries(fc)) {
            filesHtml += `<span class="file-badge">${dir}: ${count}</span>`;
        }

        // Figures gallery (thumbnails)
        let figuresHtml = '';
        if (s.figures && s.figures.length > 0) {
            for (const fig of s.figures.slice(0, 12)) {
                const name = fig.split('/').pop();
                if (name.match(/\.(png|jpg|jpeg|gif|svg)$/i)) {
                    figuresHtml += `<div class="fig-thumb" title="${fig}"><span>${name}</span></div>`;
                }
            }
        }
        if (!figuresHtml) {
            figuresHtml = '<span class="no-figures">No figures yet</span>';
        }

        this.container.innerHTML = `
            <div class="lab-header">
                <div class="lab-title">Autonomous Lab</div>
                <div class="lab-iteration">Iteration ${s.iteration}</div>
                <div class="lab-role ${roleClass}">Next: ${roleLabel}</div>
                <div class="lab-status">${s.status}</div>
            </div>
            <div class="lab-body">
                <div class="lab-section">
                    <h4>Paper Progress</h4>
                    <div class="paper-progress">${paperHtml}</div>
                </div>
                <div class="lab-section">
                    <h4>Figures</h4>
                    <div class="figures-gallery">${figuresHtml}</div>
                </div>
                <div class="lab-section">
                    <h4>Files</h4>
                    <div class="file-counts">${filesHtml}</div>
                </div>
            </div>
        `;
    }
}

// Export singleton
window.LabDashboard = new LabDashboard();
