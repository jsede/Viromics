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
    console.error('Dashboard component crashed:', error, errorInfo);
  }

  render() {
    if (this.state.hasError) {
      return (
        <div className="p-8 glass-panel rounded-2xl border border-red-200 bg-red-50/60 text-center flex flex-col items-center gap-4 my-6">
          <div className="w-12 h-12 rounded-full bg-red-100 border border-red-200 flex items-center justify-center">
            <AlertTriangle className="w-6 h-6 text-red-600" />
          </div>
          <div>
            <h3 className="text-base font-bold text-slate-900">Component Display Error</h3>
            <p className="text-xs text-slate-600 mt-1 max-w-md">
              An error occurred while rendering this tab ({this.props.componentName || 'Active View'}).
            </p>
            {this.state.error && (
              <pre className="mt-3 p-3 bg-white border border-slate-200 rounded-xl text-left text-[11px] font-mono text-red-700 max-w-lg overflow-x-auto select-text">
                {this.state.error.message || String(this.state.error)}
              </pre>
            )}
          </div>
          <button
            onClick={() => this.setState({ hasError: false, error: null })}
            className="flex items-center gap-1.5 px-4 py-2 bg-white border border-slate-300 text-xs text-slate-800 font-semibold rounded-xl hover:bg-slate-50 transition-colors"
          >
            <RefreshCw className="w-3.5 h-3.5 text-red-600" />
            Retry View
          </button>
        </div>
      );
    }
    return this.props.children;
  }
}
