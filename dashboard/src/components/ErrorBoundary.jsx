import React from 'react';
import { AlertTriangle, RefreshCw } from 'lucide-react';

export default class ErrorBoundary extends React.Component {
  constructor(props) {
    super(props);
    this.state = { hasError: false, error: null };
  }

  static getDerivedStateFromError(error) {
    return { hasError: true, error };
  }

  componentDidCatch(error, errorInfo) {
    console.error("Dashboard component crashed:", error, errorInfo);
  }

  render() {
    if (this.state.hasError) {
      return (
        <div className="p-8 glass-panel rounded-2xl border border-red-500/30 bg-red-950/20 text-center flex flex-col items-center gap-4 my-6">
          <div className="w-12 h-12 rounded-full bg-red-500/10 border border-red-500/20 flex items-center justify-center">
            <AlertTriangle className="w-6 h-6 text-brand-rose" />
          </div>
          <div>
            <h3 className="text-base font-bold text-white font-outfit">Component Display Error</h3>
            <p className="text-xs text-gray-400 mt-1 max-w-md">
              An error occurred while rendering this tab ({this.props.componentName || 'Active View'}).
            </p>
            {this.state.error && (
              <pre className="mt-3 p-3 bg-dark-bg/80 border border-dark-border rounded-xl text-left text-[11px] font-mono text-brand-rose max-w-lg overflow-x-auto select-text">
                {this.state.error.message || String(this.state.error)}
              </pre>
            )}
          </div>
          <button
            onClick={() => this.setState({ hasError: false, error: null })}
            className="flex items-center gap-1.5 px-4 py-2 bg-dark-card border border-dark-border text-xs text-gray-200 font-semibold rounded-xl hover:bg-dark-hover transition-colors"
          >
            <RefreshCw className="w-3.5 h-3.5 text-brand-teal" />
            Retry View
          </button>
        </div>
      );
    }
    return this.props.children;
  }
}
